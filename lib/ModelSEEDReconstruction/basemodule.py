from __future__ import absolute_import

import logging
import os
import copy
import json
import re
from os.path import exists
from ModelSEEDReconstruction.sdkhelper import SDKHelper

logger = logging.getLogger(__name__)

class BaseModule(SDKHelper):
    def __init__(self,name,ws_client,working_dir,config):
        self.ws_client = ws_client
        self.working_dir = working_dir
        logging.basicConfig(format='%(created)s %(levelname)s: %(message)s',
                            level=logging.INFO)
        self.obj_created = []
        self.input_objects = []
        self.method = None
        self.params = None
        self.name = name
        self.config = config
        self.initialized = False
        self.ws_id = None
        self.ws_name = None
    
    def initialize_call(self,method,params,print_params=False):
        if not self.initialized:
            self.obj_created = []
            self.input_objects = []
            self.method = method
            self.params = copy.deepcopy(params)
            self.initialized = True
            if "workspace" in params:
                self.set_ws(params["workspace"])
            if print_params:
                print(json.dumps(params,indent=4))
    
    def provenance(self):
        return [{
            'description': '',
            'input_ws_objects': self.input_objects,
            'method': self.method,
            'script_command_line': "",
            'method_params': [self.params],
            'service': self.name,
            'service_ver': self.config["version"],
            # 'time': '2015-12-15T22:58:55+0000'
        }]
    
    def set_ws(self,workspace):
        if self.ws_id == workspace or self.ws_name == workspace:
            return 
        if not isinstance(workspace, str) or re.search('^\d+$',workspace) != None:
            if isinstance(workspace, str):
                workspace = int(workspace)
            self.ws_id = workspace
            info = self.kbase_api.ws_client.get_workspace_info({"id":workspace})
            self.ws_name = info[1]
        else:
            self.ws_name = workspace
            info = self.kbase_api.ws_client.get_workspace_info({"workspace":workspace})
            self.ws_id = info[0]    