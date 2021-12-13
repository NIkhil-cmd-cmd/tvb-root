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

"""
.. moduleauthor:: Paula Popa <paula.popa@codemart.ro>
.. moduleauthor:: Bogdan Valean <bogdan.valean@codemart.ro>
"""

import os
from functools import partial

from tvb.basic.logger.builder import get_logger
from tvb.config import SIMULATOR_MODULE, SIMULATOR_CLASS
from tvb.core.adapters.abcadapter import ABCAdapter
from tvb.core.entities.load import get_class_by_name
from tvb.core.entities.model.model_operation import STATUS_ERROR, STATUS_CANCELED, STATUS_FINISHED
from tvb.core.entities.model.model_operation import STATUS_STARTED, STATUS_PENDING
from tvb.core.entities.storage import dao
from tvb.core.services.backend_clients.hpc_client import HPCJobStatus
from tvb.core.services.backend_clients.hpc_pipeline_client import HPCPipelineClient
from tvb.core.services.backend_clients.hpc_scheduler_client import HPCSchedulerClient
from tvb.core.services.exceptions import OperationException
from tvb.storage.storage_interface import StorageInterface

try:
    from pyunicore.client import Job, Transport
except ImportError:
    from tvb.basic.config.settings import HPCSettings

    HPCSettings.CAN_RUN_HPC = False


class HPCOperationService(object):
    LOGGER = get_logger(__name__)

    @staticmethod
    def _operation_error(operation, auth_token):
        operation = dao.get_operation_by_id(operation.id)
        algorithm = operation.algorithm
        adapter_instance = ABCAdapter.build_adapter(algorithm)
        if type(adapter_instance) is get_class_by_name("{}.{}".format(SIMULATOR_MODULE, SIMULATOR_CLASS)):
            simulator_gid = operation.view_model_gid
            HPCOperationService._operation_finished(operation, simulator_gid)
        else:
            op_idents = dao.get_operation_process_for_operation(operation.id, True)
            wds = []
            for op_ident in op_idents:
                job = Job(Transport(auth_token), op_ident.job_id)
                wds.append(job.working_dir)

            operation = dao.get_operation_by_id(operation.id)
            # folder = HPCSchedulerClient.storage_interface.get_project_folder(operation.project.name)
            # storage_interface = StorageInterface()
            # if storage_interface.encryption_enabled():
            #     storage_interface.inc_project_usage_count(folder)
            #     storage_interface.sync_folders(folder)

            pipeline_results_zip = HPCPipelineClient.stage_out_logs(wds, operation)

            operation.mark_complete(STATUS_ERROR)
            dao.store_entity(operation)
            HPCPipelineClient().update_db_with_results(operation, pipeline_results_zip)

    @staticmethod
    def _operation_canceled(operation):
        operation.mark_complete(STATUS_CANCELED)
        dao.store_entity(operation)

    @staticmethod
    def _operation_started(operation):
        operation.start_now()
        dao.store_entity(operation)

    @staticmethod
    def _operation_finished(operation, simulator_gid):
        op_ident = dao.get_operation_process_for_operation(operation.id)
        # TODO: Handle login
        job = Job(Transport(os.environ[HPCSchedulerClient.CSCS_LOGIN_TOKEN_ENV_KEY]),
                  op_ident.job_id)

        operation = dao.get_operation_by_id(operation.id)
        folder = HPCSchedulerClient.storage_interface.get_project_folder(operation.project.name)
        storage_interface = StorageInterface()
        if storage_interface.encryption_enabled():
            storage_interface.inc_project_usage_count(folder)
            storage_interface.sync_folders(folder)

        try:
            sim_h5_filenames, metric_op, metric_h5_filename = \
                HPCSchedulerClient.stage_out_to_operation_folder(job.working_dir, operation, simulator_gid)

            operation.mark_complete(STATUS_FINISHED)
            dao.store_entity(operation)
            HPCSchedulerClient().update_db_with_results(operation, sim_h5_filenames, metric_op, metric_h5_filename)

        except OperationException as exception:
            HPCOperationService.LOGGER.error(exception)
            HPCOperationService._operation_error(operation, os.environ[HPCSchedulerClient.CSCS_LOGIN_TOKEN_ENV_KEY])

        finally:
            if storage_interface.encryption_enabled():
                storage_interface.sync_folders(folder)
                storage_interface.set_project_inactive(operation.project)

    @staticmethod
    def _pipeline_operation_finished(operation, auth_token):
        op_ident = dao.get_operation_process_for_operation(operation.id)
        # TODO: Handle login
        job = Job(Transport(auth_token), op_ident.job_id)

        operation = dao.get_operation_by_id(operation.id)
        # folder = HPCSchedulerClient.storage_interface.get_project_folder(operation.project.name)
        # storage_interface = StorageInterface()
        # if storage_interface.encryption_enabled():
        #     storage_interface.inc_project_usage_count(folder)
        #     storage_interface.sync_folders(folder)

        try:
            pipeline_results_zip = HPCPipelineClient.stage_out_to_operation_folder(job.working_dir, operation)

            operation.mark_complete(STATUS_FINISHED)
            dao.store_entity(operation)
            HPCPipelineClient().update_db_with_results(operation, pipeline_results_zip)

        except OperationException as exception:
            HPCOperationService.LOGGER.error(exception)
            HPCOperationService._operation_error(operation, auth_token)

        # finally:
        #     if storage_interface.encryption_enabled():
        #         storage_interface.sync_folders(folder)
        #         storage_interface.set_project_inactive(operation.project)

    @staticmethod
    def handle_hpc_status_changed(operation, simulator_gid, new_status):
        # type: (Operation, str, str) -> None

        switcher = {
            STATUS_ERROR: HPCOperationService._operation_error,
            STATUS_CANCELED: HPCOperationService._operation_canceled,
            STATUS_STARTED: HPCOperationService._operation_started,
            STATUS_FINISHED: partial(HPCOperationService._operation_finished, simulator_gid=simulator_gid)
        }
        update_func = switcher.get(new_status, lambda: "Invalid operation status")
        update_func(operation)

    @staticmethod
    def stop_operation(operation, auth_token):
        operation = dao.get_operation_by_id(operation.id)
        algorithm = operation.algorithm
        adapter_instance = ABCAdapter.build_adapter(algorithm)
        if type(adapter_instance) is get_class_by_name("{}.{}".format(SIMULATOR_MODULE, SIMULATOR_CLASS)):
            simulator_gid = operation.view_model_gid
            HPCOperationService._operation_finished(operation, simulator_gid)
        else:
            HPCOperationService._pipeline_operation_finished(operation, auth_token)

    @staticmethod
    def check_operations_job(auth_token=HPCSchedulerClient.CSCS_LOGIN_TOKEN_ENV_KEY, algos=None, current_user_id=None):
        operations = dao.get_operations_for_hpc_job(algos, current_user_id)
        if operations is None or len(operations) == 0:
            return

        for operation in operations:
            HPCOperationService.LOGGER.info("Start processing operation {}".format(operation.id))
            try:
                op_ident = dao.get_operation_process_for_operation(operation.id)
                if op_ident is not None:
                    transport = Transport(auth_token)
                    job = Job(transport, op_ident.job_id)
                    job_status = job.properties['status']
                    if job.is_running():
                        if operation.status == STATUS_PENDING and job_status == HPCJobStatus.READY.value:
                            HPCOperationService._operation_started(operation)
                        HPCOperationService.LOGGER.info(
                            "CSCS job status: {} for operation {}.".format(job_status, operation.id))
                        return
                    HPCOperationService.LOGGER.info(
                        "Job for operation {} has status {}".format(operation.id, job_status))
                    if job_status == HPCJobStatus.SUCCESSFUL.value:
                        HPCOperationService.stop_operation(operation, auth_token)
                    else:
                        HPCOperationService._operation_error(operation, auth_token)
            except Exception:
                HPCOperationService.LOGGER.error(
                    "There was an error on background processing process for operation {}".format(operation.id),
                    exc_info=True)

    @staticmethod
    def check_pipeline_operations_job(auth_token=HPCSchedulerClient.CSCS_LOGIN_TOKEN_ENV_KEY, algos=None,
                                      current_user_id=None):
        operations = dao.get_operations_for_hpc_job(algos, current_user_id)
        if operations is None or len(operations) == 0:
            return

        for operation in operations:
            HPCOperationService.LOGGER.info("Start processing operation {}".format(operation.id))
            try:
                op_identifiers = dao.get_operation_process_for_operation(operation.id, True)
                if op_identifiers:
                    if len(op_identifiers) < 2:
                        HPCOperationService.LOGGER.info("For the moment {} jobs started for operation {}.".format(
                            len(op_identifiers), operation.id))
                        return
                    transport = Transport(auth_token)
                    job1 = Job(transport, op_identifiers[0].job_id)
                    job2 = Job(transport, op_identifiers[1].job_id)
                    job_status1 = job1.properties['status']
                    job_status2 = job2.properties['status']
                    if job1.is_running() or job2.is_running():
                        if operation.status == STATUS_PENDING and (
                                job_status1 == HPCJobStatus.READY.value or job_status2 == HPCJobStatus.READY.value):
                            HPCOperationService._operation_started(operation)
                        HPCOperationService.LOGGER.info(
                            "CSCS job1 status: {} for operation {}.".format(job_status1, operation.id))
                        HPCOperationService.LOGGER.info(
                            "CSCS job2 status: {} for operation {}.".format(job_status2, operation.id))
                        return

                    HPCOperationService.LOGGER.info(
                        "Job1 for operation {} has status {}".format(operation.id, job_status1))
                    HPCOperationService.LOGGER.info(
                        "Job2 for operation {} has status {}".format(operation.id, job_status2))
                    if job_status1 == HPCJobStatus.SUCCESSFUL.value or job_status1 == HPCJobStatus.SUCCESSFUL.value:
                        HPCOperationService.stop_operation(operation, auth_token)
                    else:
                        HPCOperationService._operation_error(operation, auth_token)
            except Exception:
                HPCOperationService.LOGGER.error(
                    "There was an error on background processing process for operation {}".format(operation.id),
                    exc_info=True)
