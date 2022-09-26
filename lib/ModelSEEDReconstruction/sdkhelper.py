from __future__ import absolute_import

import os
import sys
import uuid
import logging
import jinja2

logger = logging.getLogger(__name__)

class SDKHelper:
    
    @staticmethod
    def extend_sys_paths(paths):
        paths = paths.split(";")
        for path in paths:
            sys.path.append(path)
            
    @staticmethod
    def validate_args(params,required,defaults):
        #print("One:"+json.dumps(params,indent=4)+"\n\n")
        for item in required:
            if item not in params:
                raise ValueError('Required argument '+item+' is missing!')
        for key in defaults:
            if key not in params:
                params[key] = defaults[key]
        return params
    
    @staticmethod
    def create_ref(id,workspace):
        if isinstance(id, int):
            id = str(id)
        if isinstance(workspace, int):
            workspace = str(workspace)    
        return workspace+"/"+id
        
    @staticmethod
    def process_ws_ids(id_or_ref, workspace=None,no_ref=False):
        """
        IDs should always be processed through this function so we can interchangeably use
        refs, IDs, and names for workspaces and objects
        """
        objspec = {}
        if len(id_or_ref.split("/")) > 1:
            if no_ref:
                array = id_or_ref.split("/")
                workspace = array[0]
                id_or_ref = array[1]
            else:
                objspec["ref"] = id_or_ref
        if "ref" not in objspec:
            if isinstance(workspace, int):
                objspec['wsid'] = workspace
            else:
                objspec['workspace'] = workspace
            if isinstance(id_or_ref, int):
                objspec['objid'] = id_or_ref
            else:
                objspec['name'] = id_or_ref
        return objspec
    
    
    
    @staticmethod
    def build_report(context,template_file):
        # Directory this file is in
        array = template_file.split("/")
        filename = array.pop()
        template_dir = "/".join(array)
        env = jinja2.Environment(
                loader=jinja2.FileSystemLoader(template_dir),
                autoescape=jinja2.select_autoescape(['html', 'xml']))
        # Return string of html
        return env.get_template(filename).render(context)