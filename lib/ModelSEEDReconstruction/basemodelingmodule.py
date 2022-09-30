from __future__ import absolute_import

import logging
import os
import sys
import json
import cobrakbase
from ModelSEEDReconstruction.basemodule import BaseModule
from cobrakbase.core.kbasefba.fbamodel_from_cobra import CobraModelConverter
from os.path import exists

logger = logging.getLogger(__name__)

class BaseModelingModule(BaseModule):
    def __init__(self,name,ws_client,working_dir,config):
        BaseModule.__init__(self,name,ws_client,working_dir,config)
        logging.basicConfig(format='%(created)s %(levelname)s: %(message)s',
                            level=logging.INFO)
        self.kbase_api = cobrakbase.KBaseAPI(token="...")
        self.kbase_api.ws_client = ws_client
    
    def save_model(self,mdlutl,workspace=None):
        if workspace:
            self.set_ws(workspace)
        fbamodel = CobraModelConverter(mdlutl.model, mdlutl.model.genome, mdlutl.model.template).build()  # later assign core template to None
        json = fbamodel.get_data()
        mdlutl.create_kb_gapfilling_data(json,self.config["ATP_media_workspace"])
        params = {
            'id':self.ws_id,
            'objects': [{
                'data': json,
                'name': mdlutl.model.id,
                'type': "KBaseFBA.FBAModel",
                'meta': {},
                'provenance': self.provenance()
            }]
        }
        self.ws_client.save_objects(params)
        self.obj_created.append({"ref":self.create_ref(mdlutl.model.id,self.ws_name),"description":""})