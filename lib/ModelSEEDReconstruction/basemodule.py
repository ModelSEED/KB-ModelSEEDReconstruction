# -*- coding: utf-8 -*-

from __future__ import absolute_import

import os
import uuid
import logging
from cobrakbase.sdk.sdkhelper import SDKHelper
from installed_clients.KBaseReportClient import KBaseReport
from installed_clients.DataFileUtilClient import DataFileUtil
from installed_clients.WorkspaceClient import Workspace as Workspace

logger = logging.getLogger(__name__)

class SDKSupportModule:
    def __init__(self,config,callback_url = None):
        self.config = config
        self.dfu = DataFileUtil(self.callback_url)
        if callback_url == None and "SDK_CALLBACK_URL" in os.environ:
            self.callback_url = os.environ['SDK_CALLBACK_URL']
        logging.basicConfig(format='%(created)s %(levelname)s: %(message)s',
            level=logging.INFO)
        self.scratch = self.config["scratch"]
        self.report_info = None
        self.report_html = None
        self.output_type = None
        self.output_id = None
        self.ctx = None
        self.wsclient = None
        self.workspace = None
        self.objects_created = []
    
    def clear_context(self):
        self.report_info = None
        self.ctx = None
        self.output_type = None
        self.output_id = None
        self.wsclient = None
    
    def finalize_call(self,output):
        if self.report_info != None:
            output['report_name'] = self.report_info['name']
            output['report_ref'] = self.report_info['ref']
        if self.workspace != None:   
            output['workspace_name'] = self.workspace
            output['ws'] = self.workspace
        if self.output_type != None:
            output['type'] = self.output_type
            output['obj'] = self.output_id
        return output

    def initialize_call(self,ctx,workspace=None,output_type = None,output_id = None):
        self.clear_context()
        self.workspace = workspace
        self.ctx = ctx
        self.output_type = output_type
        self.output_id = output_id
        self.objects_created = []
        self.wsclient = Workspace(self.config["workspace-url"], token=self.ctx['token'])
    
    def add_created_object(self,ref,description):
        self.objects_created.append({"ref":ref,"description":description})

    def save_to_ws(self,object,type,id_or_ref,workspace = None,meta = None,hidden = False,description = ""):
        objspec = SDKHelper.process_ws_ids(id_or_ref,workspace,True)
        ws_params = {'objects':[{'type': type,'provenance': ctx["provenance"],'data': object}]}
        if "wsid" in objspec:
            ws_params["id"] = objspec["wsid"]
        else:
            ws_params["workspace"] = objspec["workspace"]   
        if "objid" in objspec:
            ws_params["objects"][0]["objid"] = objspec["objid"]
        else:
            ws_params["objects"][0]["name"] = objspec["name"] 
        if hidden
            ws_params["objects"][0]["hidden"] = 1
        if meta != None:
            ws_params["objects"][0]["meta"] = meta
        save_output = self.wsclient.save_objects(ws_params)
        self.add_created_object(SDKHelper.info_to_ref(save_output),description)
        return save_output
    
    def create_report(self,context,report_html,height=500):
        html_report_folder = os.path.join(self.scratch, 'htmlreport')
        os.makedirs(html_report_folder, exist_ok=True)
        
        with open(os.path.join(html_report_folder, 'view.html'), 'w') as f:
            f.write(self.report_html)

        report_shock_id = ""
        if self.config["save_report_to_kbase"] == "1":
            report_shock_id = self.dfu.file_to_shock({'file_path': html_report_folder,'pack': 'zip'})['shock_id']

        html_output = {
            'name' : 'view.html',
            'shock_id': report_shock_id
        }
        report_params = {
            'objects_created': self.objects_created,
            'workspace_name': self.workspace,
            'html_links': [html_output],
            'direct_html_link_index': 0,
            'html_window_height': height,
            'report_object_name': self.name + '_report_' + str(uuid.uuid4())
        }
        if self.config["save_report_to_kbase"] == "1":
            report = KBaseReport(self.callback_url, token=self.ctx['token'])
            self.report_info = report.create_extended_report(report_params)
        return self.report_html