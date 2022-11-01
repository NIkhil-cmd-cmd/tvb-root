# -*- coding: utf-8 -*-
#
#
# TheVirtualBrain-Framework Package. This package holds all Data Management, and
# Web-UI helpful to run brain-simulations. To use it, you also need do download
# TheVirtualBrain-Scientific Package (for simulators). See content of the
# documentation-folder for more details. See also http://www.thevirtualbrain.org
#
# (c) 2012-2022, Baycrest Centre for Geriatric Care ("Baycrest") and others
#
# This program is free software: you can redistribute it and/or modify it under the
# terms of the GNU General Public License as published by the Free Software Foundation,
# either version 3 of the License, or (at your option) any later version.
# This program is distributed in the hope that it will be useful, but WITHOUT ANY
# WARRANTY; without even the implied warranty of MERCHANTABILITY or FITNESS FOR A
# PARTICULAR PURPOSE.  See the GNU General Public License for more details.
# You should have received a copy of the GNU General Public License along with this
# program.  If not, see <http://www.gnu.org/licenses/>.
#
#
#   CITATION:
# When using The Virtual Brain for scientific publications, please cite it as follows:
#
#   Paula Sanz Leon, Stuart A. Knock, M. Marmaduke Woodman, Lia Domide,
#   Jochen Mersmann, Anthony R. McIntosh, Viktor Jirsa (2013)
#       The Virtual Brain: a simulator of primate brain network dynamics.
#   Frontiers in Neuroinformatics (7:10. doi: 10.3389/fninf.2013.00010)
#
#

import os
import time
from contextlib import closing
from enum import Enum
from threading import Thread, Event

from requests import HTTPError
from tvb.basic.config.settings import HPCSettings
from tvb.basic.logger.builder import get_logger
from tvb.basic.profile import TvbProfile
from tvb.core.entities.model.model_operation import OperationProcessIdentifier
from tvb.core.entities.storage import dao
from tvb.core.services.backend_clients.backend_client import BackendClient
from tvb.core.services.exceptions import OperationException
from tvb.storage.storage_interface import StorageInterface

try:
    import pyunicore.client as unicore_client
    from pyunicore.client import Job, Storage, Client
except ImportError:
    HPCSettings.CAN_RUN_HPC = False

LOGGER = get_logger(__name__)

HPC_THREADS = []


class HPCJobStatus(Enum):
    STAGINGIN = "STAGINGIN"
    READY = "READY"
    QUEUED = "QUEUED"
    STAGINGOUT = "STAGINGOUT"
    SUCCESSFUL = "SUCCESSFUL"
    FAILED = "FAILED"


def get_op_thread(op_id):
    # type: (int) -> HPCOperationThread
    op_thread = None
    for thread in HPC_THREADS:
        if thread.operation_id == op_id:
            op_thread = thread
            break
    if op_thread is not None:
        HPC_THREADS.remove(op_thread)
    return op_thread


class HPCOperationThread(Thread):
    def __init__(self, operation_id, *args, **kwargs):
        super(HPCOperationThread, self).__init__(*args, **kwargs)
        self.operation_id = operation_id
        self._stop_event = Event()

    def stop(self):
        self._stop_event.set()

    def stopped(self):
        return self._stop_event.is_set()


class HPCClientBase(BackendClient):
    """
    Base class for HPC clients: HPCPipelineClient and HPCSchedulerClient (simulation).
    """
    OUTPUT_FOLDER = 'output'
    TVB_BIN_ENV_KEY = 'TVB_BIN'
    # TODO: handle this differently in subclasses
    CSCS_LOGIN_TOKEN_ENV_KEY = 'CSCS_LOGIN_TOKEN'
    CSCS_PROJECT = 'CSCS_PROJECT'
    HOME_FOLDER_MOUNT = '/HOME_FOLDER'
    CSCS_DATA_FOLDER = 'data'
    CONTAINER_INPUT_FOLDER = '/home/tvb_user/.data'
    storage_interface = StorageInterface()

    @staticmethod
    def compute_available_disk_space(operation):
        # type: (Operation) -> int
        disk_space_per_user = TvbProfile.current.MAX_DISK_SPACE
        pending_op_disk_space = dao.compute_disk_size_for_started_ops(operation.fk_launched_by)
        user_disk_space = dao.compute_user_generated_disk_size(operation.fk_launched_by)  # From kB to Bytes
        available_space = disk_space_per_user - pending_op_disk_space - user_disk_space
        return available_space

    @staticmethod
    def _prepare_pyunicore_job(operation, job_inputs, job_script, job_config, auth_token,
                               inputs_subfolder=''):
        # use "DAINT-CSCS" -- change if another supercomputer is prepared for usage
        LOGGER.info("Prepare unicore client for operation: {}".format(operation.id))
        site_client = HPCClientBase._build_unicore_client(auth_token,
                                                          unicore_client._HBP_REGISTRY_URL,
                                                          TvbProfile.current.hpc.HPC_COMPUTE_SITE)

        LOGGER.info("Submit job for operation: {}".format(operation.id))
        job = HPCClientBase._create_job_with_pyunicore(pyunicore_client=site_client, job_description=job_config,
                                                       job_script=job_script, inputs=job_inputs,
                                                       inputs_subfolder=inputs_subfolder)
        LOGGER.info("Job url {} for operation: {}".format(job.resource_url, operation.id))
        op_identifier = OperationProcessIdentifier(operation_id=operation.id, job_id=job.resource_url)
        dao.store_entity(op_identifier)

        mount_point = job.working_dir.properties[HPCSettings.JOB_MOUNT_POINT_KEY]
        LOGGER.info("Job mount point: {}".format(mount_point))
        return job

    @staticmethod
    def _build_unicore_client(auth_token, registry_url, supercomputer):
        # type: (str, str, str) -> Client
        transport = unicore_client.Transport(auth_token)
        registry = unicore_client.Registry(transport, registry_url)
        return registry.site(supercomputer)

    @staticmethod
    def _poll_job(job):
        '''wait until job completes'''
        while job.is_running():
            time.sleep(5.1)

    @staticmethod
    def _create_job_with_pyunicore(pyunicore_client, job_description, job_script, inputs, inputs_subfolder=''):
        # type: (Client, {}, str, list, str) -> Job
        """
        Submit and start a batch job on the site, optionally uploading input data files.
        We took this code from the pyunicore Client.new_job method in order to use our own upload method
        :return: job
        """

        if len(inputs) > 0 or job_description.get('haveClientStageIn') is True:
            job_description['haveClientStageIn'] = "true"

        with closing(
                pyunicore_client.transport.post(url=pyunicore_client.links['jobs'], json=job_description)) as resp:
            job_url = resp.headers['Location']

        job = Job(pyunicore_client.transport, job_url)

        working_dir = job.working_dir
        if job_script is not None:
            HPCClientBase._upload_file_with_pyunicore(working_dir, job_script, None)
        for input_file in inputs:
            HPCClientBase._upload_file_with_pyunicore(working_dir, input_file, inputs_subfolder)
        return job

    @staticmethod
    def _upload_file_with_pyunicore(working_dir, input_name, subfolder='', destination=None):
        # type: (Storage, str, object, str) -> None
        """
        Upload file to the HPC working dir.
        We took this upload code from pyunicore Storage.upload method and modified it because in the original code the
        upload URL is generated using the os.path.join method. The result is an invalid URL for windows os.
        """
        if destination is None:
            destination = os.path.basename(input_name)
        if subfolder == '':
            subfolder = HPCClientBase.CSCS_DATA_FOLDER

        if subfolder:
            url = "{}/{}/{}/{}".format(working_dir.resource_url, "files", subfolder, destination)
        else:
            url = "{}/{}/{}".format(working_dir.resource_url, "files", destination)

        headers = {'Content-Type': 'application/octet-stream'}
        with open(input_name, 'rb') as fd:
            working_dir.transport.put(
                url=url,
                headers=headers,
                data=fd)

    @staticmethod
    def _listdir(working_dir, base='/'):
        # type: (Storage, str) -> dict
        """
        We took this code from pyunicore Storage.listdir method and extended it to use a subdirectory.
        Looking at the method signature, it should have had this behavior, but the 'base' argument is not used later
        inside the method code.
        Probably will be fixed soon in their API, so we could delete this.
        :return: dict of {str: PathFile} objects
        """
        ret = {}
        try:
            for path, meta in working_dir.contents(base)['content'].items():
                path_url = working_dir.path_urls['files'] + path
                path = path[1:]  # strip leading '/'
                if meta['isDirectory']:
                    ret[path] = unicore_client.PathDir(working_dir, path_url, path)
                else:
                    ret[path] = unicore_client.PathFile(working_dir, path_url, path)
            return ret
        except HTTPError as http_error:
            if http_error.response.status_code == 404:
                raise OperationException("Folder {} is not present on HPC storage.".format(base))
            raise http_error

    @staticmethod
    def _stage_out_outputs(encrypted_dir_path, output_list):
        # type: (str, dict) -> list
        if not os.path.isdir(encrypted_dir_path):
            os.makedirs(encrypted_dir_path)

        encrypted_files = list()
        for output_filename, output_filepath in output_list.items():
            if type(output_filepath) is not unicore_client.PathFile:
                LOGGER.info("Object {} is not a file.")
                continue
            filename = os.path.join(encrypted_dir_path, os.path.basename(output_filename))
            output_filepath.download(filename)
            encrypted_files.append(filename)
        return encrypted_files
