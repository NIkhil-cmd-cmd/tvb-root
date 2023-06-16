import requests
import json
import pooch
from pathlib import Path
import logging
from zipfile import ZipFile

from .base import ZenodoDataset
from .zenodo import Zenodo, Record, BASE_URL

class TVB_Data(ZenodoDataset):

    CONCEPTID = "3417206"
    
    def __init__(self, version= "2.7"):
        """
        Constructor for TVB_Data class 

        parameters
        -----------

        version: str
              - Version number of the dataset, Default value is 2.7

        """
        super().__init__(version)
        try:
            self.recid = self.read_cached_response()[version]['conceptrecid']
        except:
            logging.warning("Data not found in cached response, updating the cached responses")
            self.recid = Zenodo().get_versions_info(self.CONCEPTID)[version]            
            self.update_cached_response() 
        
        self.rec = Record(self.read_cached_response()[self.version])
        
        logging.info(f"instantiated TVB_Data class with version {version}")
    
    def download(self):
        """
        Downloads the dataset to the cached location, skips download is file already present at the path.
        """
        self.rec.download()

    def fetch_data(self, file_name=None, extract_dir=None):        
        """
        Fetches the data 

        parameters:
        -----------
        file_name: str
                - Name of the file from the downloaded zip file to fetch. If `None`, extracts whole archive. Default is `None` 
        extract_dir: str
                - Path where you want to extract the archive, if `None` extracts the archive to current working directory. Default is `None` 


        returns: Pathlib.Path
            path of the file which was extracted
        """

        #TODO: errrors when absolute path given. 
        try:
            file_path = self.rec.file_loc['tvb_data.zip']
        except:
            self.download()
            file_path = self.rec.file_loc['tvb_data.zip']

 
        if file_name == None:
            ZipFile(file_path).extractall(path=extract_dir)
            if extract_dir==None:
                return Path.cwd()
            return Path.cwd()/ Path(extract_dir)

        with ZipFile(file_path) as zf:
            file_names_in_zip = zf.namelist()
        zf.close()

        file_names_in_zip = {str(Path(i).name): i for i in file_names_in_zip}
        ZipFile(file_path).extract(file_names_in_zip[file_name])
        return Path.cwd() / file_names_in_zip[file_name]


    def update_cached_response(self):
        """
        gets responses from zenodo server and saves them to cache file. 
        """
        
        file_dir = pooch.os_cache("pooch")/ "tvb_cached_responses.txt"
        
        responses = {}

        url = f"{BASE_URL}records?q=conceptrecid:{self.CONCEPTID}&all_versions=true"

        for hit in requests.get(url).json()['hits']['hits']:
            version = hit['metadata']['version']
            response = hit 

            responses[version] = response

        with open(file_dir, "w") as fp:
            json.dump(responses, fp)
        fp.close()

        return 

    def read_cached_response(self):
        """
        reads responses from the cache file.

        """
        
        file_dir = pooch.os_cache("pooch") / "tvb_cached_responses.txt"


        with open(file_dir) as fp:
            responses = json.load(fp)

        fp.close()


        responses = dict(responses)
        return responses
