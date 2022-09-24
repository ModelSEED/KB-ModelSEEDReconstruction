# -*- coding: utf-8 -*-
import os
import time
import unittest
from configparser import ConfigParser

from ModelSEEDReconstruction.ModelSEEDReconstructionImpl import ModelSEEDReconstruction
from ModelSEEDReconstruction.ModelSEEDReconstructionServer import MethodContext
from ModelSEEDReconstruction.authclient import KBaseAuth as _KBaseAuth

from installed_clients.WorkspaceClient import Workspace


class ModelSEEDReconstructionTest(unittest.TestCase):

    @classmethod
    def setUpClass(cls):
        token = os.environ.get('KB_AUTH_TOKEN', None)
        config_file = os.environ.get('KB_DEPLOYMENT_CONFIG', None)
        cls.cfg = {}
        config = ConfigParser()
        config.read(config_file)
        for nameval in config.items('ModelSEEDReconstruction'):
            cls.cfg[nameval[0]] = nameval[1]
        # Getting username from Auth profile for token
        authServiceUrl = cls.cfg['auth-service-url']
        auth_client = _KBaseAuth(authServiceUrl)
        user_id = auth_client.get_user(token)
        # WARNING: don't call any logging methods on the context object,
        # it'll result in a NoneType error
        cls.ctx = MethodContext(None)
        cls.ctx.update({'token': token,
                        'user_id': user_id,
                        'provenance': [
                            {'service': 'ModelSEEDReconstruction',
                             'method': 'please_never_use_it_in_production',
                             'method_params': []
                             }],
                        'authenticated': 1})
        cls.wsURL = cls.cfg['workspace-url']
        cls.wsClient = Workspace(cls.wsURL)
        cls.serviceImpl = ModelSEEDReconstruction(cls.cfg)
        cls.scratch = cls.cfg['scratch']
        cls.callback_url = os.environ['SDK_CALLBACK_URL']
        suffix = int(time.time() * 1000)
        cls.wsName = "test_ContigFilter_" + str(suffix)
        ret = cls.wsClient.create_workspace({'workspace': cls.wsName})  # noqa

    @classmethod
    def tearDownClass(cls):
        if hasattr(cls, 'wsName'):
            cls.wsClient.delete_workspace({'workspace': cls.wsName})
            print('Test workspace was deleted')

    # NOTE: According to Python unittest naming rules test method names should start from 'test'. # noqa
    def test_your_method(self):
        #All input data is stored in the following public narrative in appdev:
        #https://appdev.kbase.us/narrative/68303
        #Set output WS - initialy set to fixed chenry workspace so we can view the output in narrative:
        #https://appdev.kbase.us/narrative/68304
        #If you are not chenry, either uncomment the following line, OR ask chenry for write permission to 68304
        #output_ws = self.wsName()
        output_ws = 68304
        ret = self.serviceImpl.build_metabolic_models(self.ctx,{
            "workspace":output_ws,
            "genome_refs":["68303/Escherichia_coli_K-12_MG1655.RAST"],
            "run_gapfilling":False,
            "atp_safe":True,
            "forced_atp_list":[],
            "gapfilling_media_list":["KBaseMedia/Carbon-D-Glucose"],
            "suffix":".mdl",
            "core_template":"auto",
            "gs_template":"auto",
            "gs_template_ref":None,
            "template_reactions_only":True,
            "output_core_models":True
        })
        
        ret = self.serviceImpl.build_metabolic_models(self.ctx,{
            "workspace":output_ws,
            "genome_refs":["68303/Escherichia_coli_K-12_MG1655.RAST"],
            "run_gapfilling":True,
            "atp_safe":True,
            "forced_atp_list":[],
            "gapfilling_media_list":["KBaseMedia/Carbon-D-Glucose"],
            "suffix":".mdl2",
            "core_template":"auto",
            "gs_template":"auto",
            "gs_template_ref":None,
            "template_reactions_only":True,
            "output_core_models":True
        })
        # Check returned data with
        # self.assertEqual(ret[...], ...) or other unittest methods