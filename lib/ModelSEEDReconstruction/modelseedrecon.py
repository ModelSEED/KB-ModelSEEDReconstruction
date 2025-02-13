from __future__ import absolute_import

import os
import sys
import uuid
import logging
import json
import jinja2
import pandas as pd
from optlang.symbolics import Zero, add
import cobra
from cobrakbase.core.kbasefba import FBAModel
from modelseedpy import AnnotationOntology, MSPackageManager,MSGenome, MSMedia, MSModelUtil, MSBuilder, MSGapfill, FBAHelper, MSGrowthPhenotypes, MSModelUtil, MSATPCorrection,MSModelReport
from modelseedpy.helpers import get_template
from modelseedpy.core.msgenomeclassifier import MSGenomeClassifier
from modelseedpy.core.mstemplate import MSTemplateBuilder
from modelseedpy.core.msgenome import normalize_role
from kbbasemodules.basemodelingmodule import BaseModelingModule
from modelseedpy.community.mscommunity import MSCommunity

logger = logging.getLogger(__name__)

class ModelSEEDRecon(BaseModelingModule):
    def __init__(self,config,module_dir="/kb/module",working_dir=None,token=None,clients={},callback=None):
        BaseModelingModule.__init__(self,"ModelSEEDReconstruction",config,module_dir=module_dir,working_dir=working_dir,token=token,clients=clients,callback=callback)
        self.core_template = None
        self.util = None
        self.gs_template = None
        self.version = "0.1.1.msr"
        self.module_dir = module_dir
        self.native_ontology = False
        logging.basicConfig(format='%(created)s %(levelname)s: %(message)s',
                            level=logging.INFO)
        
    def build_metabolic_models(self,params):
        default_media = "KBaseMedia/AuxoMedia"
        self.initialize_call("build_metabolic_models",params,True)
        if "genomes" in params:
            params["genome_refs"] = params["genomes"]["genome_refs"]
        if "reconstruction_parameters" in params:
            params["suffix"] = params["reconstruction_parameters"]["suffix"]
            #params["output_core_models"] = params["reconstruction_parameters"]["output_core_models"]
        if "gapfilling_parameters" in params:
            params["run_gapfilling"] = params["gapfilling_parameters"]["run_gapfilling"]
            params["gapfilling_media_list"] = params["gapfilling_parameters"]["gapfilling_media_list"]
        if "atp_safe_parameters" in params:
            params["atp_safe"] = params["atp_safe_parameters"]["atp_safe"]
            params["forced_atp_list"] = params["atp_safe_parameters"]["forced_atp_list"]
            params["automated_atp_evaluation"] = params["atp_safe_parameters"]["automated_atp_evaluation"]
        if "template_parameters" in params:
            params["gs_template"] = params["template_parameters"]["gs_template"]
            params["gs_template_ref"] = params["template_parameters"]["gs_template_ref"]
            params["core_template_ref"] = params["template_parameters"]["core_template_ref"]
        self.validate_args(params,["workspace"],{
            "genome_refs":[],
            "run_gapfilling":False,
            "atp_safe":True,
            "forced_atp_list":[],
            "gapfilling_media_list":None,
            "suffix":".mdl",
            "core_template":"auto",
            "gs_template":"auto",
            "gs_template_ref":None,
            "core_template_ref":None,
            "template_reactions_only":True,
            "output_core_models":False,
            "automated_atp_evaluation":True,
            "atp_medias":[],
            "load_default_medias":True,
            "max_gapfilling":10,
            "gapfilling_delta":0,
            "return_model_objects":False,
            "return_data":False,
            "save_report_to_kbase":True,
            "change_to_complete":False,
            "gapfilling_mode":"Sequential",
            "base_media":None,
            "compound_list":None,
            "base_media_target_element":"C",
            "expression_refs":None,
            "extend_model_with_ontology":False,
            "ontology_events":None,
            "save_models_to_kbase":True,
            "save_gapfilling_fba_to_kbase":True,
            "annotation_priority":[],
            "merge_annotations":False
        })
        if params["change_to_complete"]:
            default_media = "KBaseMedia/Complete"
        params["genome_refs"] = self.process_genome_list(params["genome_refs"],params["workspace"])
        #Processing media
        params["gapfilling_media_objs"] = self.process_media_list(params["gapfilling_media_list"],default_media,params["workspace"])
        #Preloading core and preselected template - note the default GS template is None because this signals for the classifier to be used to select the template
        self.gs_template = None
        if params["gs_template_ref"] != None and not isinstance(params["gs_template_ref"],str):
            #Allows user to directly set the GS template object by passing the object in gs_template_ref
            self.gs_template = params["gs_template_ref"]
        if params["gs_template_ref"]:
            #Allows user to set the workspace ID of the GS template to be used
            self.gs_template = self.get_template(params["gs_template_ref"],None)
        if params["core_template_ref"] != None and not isinstance(params["core_template_ref"],str):
            #Allows user to directly set the core template object by passing the object in core_template_ref
            self.core_template = params["core_template_ref"]
        elif params["core_template_ref"]:
            #Allows user to set the workspace ID of the core template to be used
            self.core_template = self.get_template(params["core_template_ref"],None)
        else:
            #Setting the default core template
            self.core_template = self.get_template(self.templates["core"],None)  
        #Initializing classifier
        genome_classifier = self.get_classifier()
        #Initializing output data tables
        result_table = pd.DataFrame({})
        default_output = {"Model":None,"Genome":None,"Genes":None,"Class":None,
                          "Model genes":None,"Reactions":None,
                          "Core GF":None,"GS GF":None,"Growth":None,"Comments":[]}
        #Retrieving genomes and building models one by one
        mdllist = []
        for i,gen_ref in enumerate(params["genome_refs"]):
            template_type = params["gs_template"]
            #Getting RAST annotated genome, which will be reannotated as needed
            genome = self.get_msgenome_from_ontology(gen_ref,native_python_api=self.native_ontology,output_ws=params["workspace"])
            #Initializing output row
            current_output = default_output.copy()
            current_output["Comments"] = []
            gid = genome.id
            current_output["Model"] = gid+params["suffix"]+'<br><a href="'+gid+params["suffix"]+'-recon.html" target="_blank">(see reconstruction report)</a><br><a href="'+gid+params["suffix"]+'-full.html" target="_blank">(see full view)</a>'
            current_output["Genome"] = genome.annoont.info[10]["Name"]
            current_output["Genes"] = genome.annoont.info[10]["Number of Protein Encoding Genes"]
            #Pulling annotation priority
            current_output["Comments"].append("Other annotation priorities not supported by this app yet. Using RAST.")
            if template_type == "auto":
                current_output["Class"] = genome_classifier.classify(genome)
                print(current_output["Class"])
                if current_output["Class"] == "P":
                    current_output["Class"] = "Gram Positive"
                    template_type = "gp"
                elif current_output["Class"] == "N" or current_output["Class"] == "--":
                    current_output["Class"] = "Gram Negative"
                    template_type = "gn"
                elif current_output["Class"] == "A":
                    current_output["Class"] = "Archaea"
                    template_type = "ar"
                elif current_output["Class"] == "C":
                    current_output["Class"] = "Cyanobacteria"
                    template_type = "cyano"
                    current_output["Comments"].append("Cyanobacteria not yet supported. Skipping genome.")
                    result_table = result_table.append(current_output, ignore_index = True)
                    continue
                else:
                    current_output["Comments"].append("Unrecognized genome class "+current_output["Class"]+". Skipping genome.")
                    result_table = result_table.append(current_output, ignore_index = True)
                    continue
            if not self.gs_template:
                self.gs_template = self.get_template(self.templates[template_type],None)
            #Building model            
            base_model = FBAModel({'id':gid+params["suffix"], 'name':genome.scientific_name})
            builder = MSBuilder(genome, self.gs_template)
            mdl = builder.build(base_model, '0', False, False)            
            mdl.genome = genome
            mdl.template = self.gs_template
            mdl.core_template_ref = str(self.core_template.info)
            mdl.genome_ref = self.wsinfo_to_ref(genome.annoont.info)
            mdl.template_ref = str(self.gs_template.info)
            current_output["Core GF"] = "NA" 
            mdlutl = MSModelUtil.get(mdl)
            if params["extend_model_with_ontology"] or len(params["annotation_priority"]) > 0:
                #Removing reactions from model that are not in the core
                remove_list = []
                for rxn in mdl.reactions:
                    if rxn.id not in self.core_template.reactions and rxn.id[0:3] != "bio":
                        remove_list.append(rxn.id)
                mdl.remove_reactions(remove_list)
                #Now extending model with selected ontology priorities
                if params["ontology_events"] == None and len(params["annotation_priority"]) > 0:
                    params["ontology_events"] = genome.annoont.get_events_from_priority_list(params["annotation_priority"])
                self.extend_model_with_other_ontologies(mdlutl,genome.annoont,builder,prioritized_event_list=params["ontology_events"],merge_all=params["merge_annotations"])
            mdlutl.save_model("base_model.json")
            genome_objs = {mdlutl:genome}
            expression_objs = None
            if params["expression_refs"]:
                expression_objs = self.get_expression_objs(params["expression_refs"],genome_objs)
            if params["atp_safe"]:
                atpcorrection = MSATPCorrection(mdlutl,self.core_template,params["atp_medias"],load_default_medias=params["load_default_medias"],max_gapfilling=params["max_gapfilling"],gapfilling_delta=params["gapfilling_delta"],forced_media=params["forced_atp_list"],default_media_path=self.module_dir+"/data/atp_medias.tsv")
                tests = atpcorrection.run_atp_correction()
                current_output["Core GF"] = len(atpcorrection.cumulative_core_gapfilling)
            #Setting the model ID so the model is saved with the correct name in KBase
            mdlutl.get_attributes()["class"] = current_output["Class"]
            mdlutl.wsid = gid+params["suffix"]
            #Running gapfilling
            current_output["GS GF"] = "NA"
            if params["run_gapfilling"]:
                self.gapfill_metabolic_models({
                    "media_objs":params["gapfilling_media_objs"],#
                    "model_objs":[mdlutl],#
                    "genome_objs":genome_objs,#
                    "expression_objs":expression_objs,#
                    "atp_safe":params["atp_safe"],#
                    "workspace":params["workspace"],#
                    "suffix":"",#
                    "default_objective":"bio1",#
                    "output_data":{mdlutl:current_output},#
                    "forced_atp_list":params["forced_atp_list"],
                    "templates":[self.gs_template],
                    "internal_call":True,
                    "gapfilling_mode":params["gapfilling_mode"],
                    "base_media":params["base_media"],
                    "compound_list":params["compound_list"],
                    "base_media_target_element":params["base_media_target_element"],
                    "save_models_to_kbase":params["save_models_to_kbase"],
                    "save_gapfilling_fba_to_kbase":params["save_gapfilling_fba_to_kbase"]
                })
            else:
                if params["save_models_to_kbase"]:
                    self.save_model(mdlutl,params["workspace"],None)
                mdlutl.model.objective = "bio1"
                mdlutl.pkgmgr.getpkg("KBaseMediaPkg").build_package(None)
                current_output["Growth"] = "Complete:"+str(mdlutl.model.slim_optimize())
            current_output["Reactions"] = mdlutl.nonexchange_reaction_count()
            current_output["Model genes"] = len(mdlutl.model.genes)
            #if params["output_core_models"]:
            #TODO: Remove noncore reactions and change biomass and change model ID and then resave
            #Filling in model output
            result_table = pd.concat([result_table, pd.DataFrame([current_output])], ignore_index=True)
            mdllist.append(mdlutl)
        output = {}
        self.build_dataframe_report(result_table,mdllist)
        if params["save_report_to_kbase"]:
            output = self.save_report_to_kbase()
        if params["return_data"]:
            output["data"] = result_table.to_json()
        if params["return_model_objects"]:
            output["model_objs"] = mdllist
        return output
    
    def gapfill_metabolic_models(self,params):
        default_media = "KBaseMedia/AuxoMedia"
        self.initialize_call("gapfill_metabolic_models",params,True)
        self.validate_args(params,["workspace"],{
            "media_list":None,
            "media_objs":None,
            "genome_objs":None,
            "expression_refs":None,
            "expression_objs":None,
            "model_list":None,
            "model_objectives":[],
            "model_objs":[],
            "atp_safe":True,
            "suffix":".gf",
            "forced_atp_list":[],
            "templates":None,
            "core_template_ref":None,
            "source_models":[],
            "limit_medias":[],
            "limit_objectives":[],
            "limit_thresholds":[],
            "is_max_limits":[],
            "minimum_objective":0.01,
            "reaction_exlusion_list":[],
            "default_objective":"bio1",
            "kbmodel_hash":{},
            "output_data":None,
            "internal_call":False,
            "atp_medias":[],
            "load_default_atp_medias":True,
            "max_atp_gapfilling":0,
            "gapfilling_delta":0,
            "return_model_objects":False,
            "return_data":False,
            "save_report_to_kbase":True,
            "change_to_complete":False,
            "gapfilling_mode":"Sequential",
            "base_media":None,
            "compound_list":None,
            "base_media_target_element":"C",
            "save_models_to_kbase":True,
            "save_gapfilling_fba_to_kbase":True
        })
        base_comments = []
        if params["change_to_complete"]:
            base_comments.append("Changing default to complete.")
            default_media = "KBaseMedia/Complete"
        result_table = pd.DataFrame({})
        default_output = {"Model":None,"Genome":None,"Genes":None,"Class":None,
            "Model genes":None,"Reactions":None,"Core GF":None,"GS GF":None,"Growth":None,"Comments":[]}
        #Retrieving models if not provided already
        if "model_objs" not in params or len(params["model_objs"]) == 0:
            params["model_objs"] = []
            for mdl_ref in params["model_list"]:
                params["model_objs"].append(self.get_model(mdl_ref))
        #Retrieving genomes if not provided already
        if not params["genome_objs"]:
            params["genome_objs"] = {}
            for mdl in params["model_objs"]:
                params["genome_objs"][mdl] = self.get_msgenome_from_ontology(mdl.model.genome_ref,native_python_api=self.native_ontology,output_ws=params["workspace"])
        #Retrieving expression data if not provided already
        if not params["expression_objs"] and params["expression_refs"]:
            params["expression_objs"] = self.get_expression_objs(params["expression_refs"],params["genome_objs"])
        #Processing media
        if not params["media_objs"]:
            params["media_objs"] = self.process_media_list(params["media_list"],default_media,params["workspace"])
        #Processing compound list
        if params["compound_list"]:
            if not params["base_media"]:
                base_comments.append("No base media provided. Ignoring compound list.")
            else:
                for cpd in params["compound_list"]:
                    newmedia = MSMedia.from_dict({cpd:100})
                    newmedia.merge(params["base_media"])
                    params["media_objs"].append(newmedia)
        #Compiling additional tests
        additional_tests = []
        for i,limit_media in enumerate(params["limit_medias"]):
            additional_tests.append({
                "objective":params["limit_objectives"][i],
                "media":self.get_media(limit_media,None),
                "is_max_threshold":params["is_max_limits"][i],
                "threshold":params["limit_thresholds"][i]
            })
        #Getting core template
        if params["core_template_ref"]:
            self.core_template = self.get_template(params["core_template_ref"],None)
        else:
            self.core_template = self.get_template(self.templates["core"],None)
        #Iterating over each model and running gapfilling
        for i,mdlutl in enumerate(params["model_objs"]):
            current_output = default_output.copy()
            current_output["Comments"] = base_comments
            current_output["Model"] = mdlutl.wsid+params["suffix"]+'<br><a href="'+mdlutl.wsid+params["suffix"]+'-recon.html" target="_blank">(see reconstruction report)</a><br><a href="'+mdlutl.wsid+params["suffix"]+'-full.html" target="_blank">(see full view)</a>'
            if params["output_data"] and mdlutl in params["output_data"]:
                current_output = params["output_data"][mdlutl]
            #Setting the objective
            if i < len(params["model_objectives"]):
                if not params["model_objectives"][i]:
                    params["model_objectives"][i] = params["default_objective"]
            else:
                params["model_objectives"].append(params["default_objective"])
            #Computing tests for ATP safe gapfilling
            if params["atp_safe"]:
                tests = mdlutl.get_atp_tests(core_template=self.core_template,atp_media_filename=self.module_dir+"/data/atp_medias.tsv",recompute=False)
                additional_tests.extend(tests)
            #Creating gapfilling object and configuring solver
            #mdlutl.model.solver = config["solver"]
            if not params["templates"]:
                params["templates"] = [self.get_template(mdlutl.model.template_ref)]
            msgapfill = MSGapfill(
                mdlutl,
                params["templates"],
                params["source_models"],
                additional_tests,
                blacklist=params["reaction_exlusion_list"],
                default_target=params["model_objectives"][i],
                minimum_obj=params["minimum_objective"],
                base_media=params["base_media"],
                base_media_target_element=params["base_media_target_element"]
            )
            #Setting reaction scores from genome
            msgapfill.reaction_scores = params["genome_objs"][mdlutl].annoont.get_reaction_gene_hash(feature_type="gene")
            if self.util:
                self.util.save("original_scores",msgapfill.reaction_scores)
            if params["expression_objs"] and mdlutl in params["expression_objs"] and mdlutl in params["genome_objs"]:
                expression_scores = msgapfill.compute_reaction_weights_from_expression_data(params["expression_objs"][mdlutl],params["genome_objs"][mdlutl].annoont)
                for rxn_id in msgapfill.reaction_scores:
                    for gene in msgapfill.reaction_scores[rxn_id]:
                        if gene in expression_scores:
                            msgapfill.reaction_scores[rxn_id][gene]["probability"] = expression_scores[gene]+0.5
                if self.util:
                    self.util.save("expression_scores",msgapfill.reaction_scores)   
            #Running gapfilling in all conditions
            mdlutl.gfutl.cumulative_gapfilling = []
            growth_array = []
            solutions = msgapfill.run_multi_gapfill(
                params["media_objs"],
                target=params["model_objectives"][i],
                default_minimum_objective=params["minimum_objective"],
                binary_check=False,
                prefilter=True,
                check_for_growth=True,
                gapfilling_mode=params["gapfilling_mode"],
                run_sensitivity_analysis=True,
                integrate_solutions=True
            )
            output_solution = None
            output_solution_media = None
            for media in params["media_objs"]:
                if media in solutions and "growth" in solutions[media]:
                    growth_array.append(media.id+":"+str(solutions[media]["growth"]))
                    if solutions[media]["growth"] > 0 and output_solution == None:
                        mdlutl.pkgmgr.getpkg("KBaseMediaPkg").build_package(media)
                        mdlutl.pkgmgr.getpkg("ElementUptakePkg").build_package({"C": 60})
                        output_solution = cobra.flux_analysis.pfba(mdlutl.model)
                        output_solution_media = media
            solution_rxn_types = ["new","reversed"]
            if output_solution and output_solution_media in solutions:
                gfsolution = solutions[output_solution_media]
                for rxn_type in solution_rxn_types:
                    for rxn_id in gfsolution[rxn_type]:
                        if gfsolution[rxn_type][rxn_id] == ">":
                            output_solution.fluxes[rxn_id] = 1000
                        else:
                            output_solution.fluxes[rxn_id] = -1000
            current_output["Growth"] = "<br>".join(growth_array)
            current_output["GS GF"] = len(mdlutl.gfutl.cumulative_gapfilling)
            current_output["Reactions"] = mdlutl.nonexchange_reaction_count()
            current_output["Model genes"] = len(mdlutl.model.genes)
            #Saving completely gapfilled model
            if params["save_models_to_kbase"]:
                self.save_model(mdlutl,params["workspace"],None,params["suffix"])
            if params["save_gapfilling_fba_to_kbase"] and output_solution:
                self.save_solution_as_fba(output_solution,mdlutl,output_solution_media,mdlutl.wsid+".fba",workspace=params["workspace"],fbamodel_ref=str(params["workspace"])+"/"+mdlutl.wsid)
            if not params["internal_call"]:
                result_table = result_table.append(current_output, ignore_index = True) 
        output = {}
        if not params["internal_call"]:
            self.build_dataframe_report(result_table,params["model_objs"])
            if params["save_report_to_kbase"]:
                output = self.save_report_to_kbase()
            if params["return_data"]:
                output["data"] = result_table.to_json()
            if params["return_model_objects"]:
                output["model_objs"] = params["model_objs"]
        return output

    def run_community_fba(self,params):
        self.initialize_call("gapfill_metabolic_models",params,True)
        self.validate_args(params,["workspace","fbamodel_id","fba_output_id"],{
            "media_id":"KBaseMedia/AuxoMedia",
            "target_reaction":"bio1",
            "expression_ref":None,
            "expression_condition":None,
            "feature_ko_list":[],
		    "reaction_ko_list":[],
            "prfba":True,
            "objective_fraction":0.9,
            "expression_coef":-10,
            "min_probability":0.05,
            "kinetics_coef":400,
            "exchange_coef":1,
            "max_c_uptake":None,
            "max_n_uptake":None,
            "max_p_uptake":None,
            "max_s_uptake":None,
            "max_o_uptake":None,
            "predict_community_composition":False,
            "return_fba_object":False,
            "return_data":False,
            "save_report_to_kbase":True,
            "save_fba_to_kbase":True
        })
        #Getting model
        mdlutl = self.get_model(params["fbamodel_id"])
        #Setting objective
        mdlutl.model.objective = params["target_reaction"]
        #Setting media
        media = self.get_media(params["media_id"],None)
        mdlutl.pkgmgr.getpkg("KBaseMediaPkg").build_package(media)
        #Setting elemental uptake constraint
        exchange_string = "None"
        elemental_hash = {}
        if params["max_c_uptake"]:
            elemental_hash["C"] = params["max_c_uptake"]
            exchange_string += "C: "+str(params["max_c_uptake"])+"; "
        if params["max_n_uptake"]:
            elemental_hash["N"] = params["max_n_uptake"]
            exchange_string += "N: "+str(params["max_n_uptake"])+"; "
        if params["max_p_uptake"]:
            elemental_hash["P"] = params["max_p_uptake"]
            exchange_string += "P: "+str(params["max_p_uptake"])+"; "
        if params["max_s_uptake"]:
            elemental_hash["S"] = params["max_s_uptake"]
            exchange_string += "S: "+str(params["max_s_uptake"])+"; "
        if params["max_o_uptake"]:
            elemental_hash["O"] = params["max_o_uptake"]
            exchange_string += "O: "+str(params["max_o_uptake"])+"; "
        if len(elemental_hash) > 0:
            mdlutl.pkgmgr.getpkg("ElementUptakePkg").build_package(elemental_hash)
        #Adding commkinetic constraints
        species_list = []
        abundance_hash = {}
        comm_model = MSCommunity(
            model=mdlutl,
            names=species_list
        )
        comm_model.set_abundance(abundance_hash)
        mdlutl.pkgmgr.getpkg("CommKineticPkg").build_package(float(params["kinetics_coef"]), comm_model)
        output = {}
        #Optimizing and constraining objective to fraction of optimum
        output["max_growth"] = mdlutl.model.slim_optimize()
        if "C" in elemental_hash:
            output["carbon_uptake"] = mdlutl.getpkg("ElementUptakePkg").variables["elements"]["C"].primal
        if str(output["max_growth"]) != "nan":
            return output
        mdlutl.pkgmgr.getpkg("ObjConstPkg").build_package(output["max_growth"] * float(params["objective_fraction"]), None)
        #Creating minimal probability objective
        coef = {}
        for rxn in comm_model.model.reactions:
            if "rxn" == rxn.id[0:3]:
                currrxn = comm_model.model.reactions.get_by_id(rxn.id)
                if bool(params["prfba"]):
                    coef.update(
                        {
                            currrxn.forward_variable: max(
                                float(params["min_probability"]), (1 - float(rxn.probability))
                            )
                        }
                    )
                    coef.update(
                        {
                            currrxn.reverse_variable: max(
                                float(params["min_probability"]), (1 - float(rxn.probability))
                            )
                        }
                    )
                else:
                    coef.update({currrxn.forward_variable: 1})
                    coef.update({currrxn.reverse_variable: 1})
            elif "EX_" == rxn.id[0:3]:
                currrxn = comm_model.model.reactions.get_by_id(rxn.id)
                coef.update({currrxn.forward_variable: float(params["exchange_coef"])})
                coef.update({currrxn.reverse_variable: float(params["exchange_coef"])})
        #Adding expression data to minimum probability objective
        """
        for clade in feature_entries:
            total = 0
            for ftr in feature_entries[clade]:
                total += feature_entries[clade][ftr][condition]
            for ftr in feature_entries[clade]:
                feature_entries[clade][ftr][condition] = feature_entries[clade][ftr][condition]/total
        for rxn in mdlutl.model.reactions:
            highest_exp = 0
            for gene in rxn.genes:
                array = gene.id.split("_")
                array.pop()
                clade = "_".join(array)
                if gene.id in feature_entries[clade] and condition in feature_entries[clade][gene.id]:
                    if feature_entries[clade][gene.id][condition] > highest_exp:
                        highest_exp = feature_entries[clade][gene.id][condition]
            if highest_exp > min_expression: 
                coef.update(
                    {
                        rxn.forward_variable: float(params["expression_coef"]) * highest_exp
                    }
                )
                coef.update(
                    {
                        rxn.reverse_variable: float(params["expression_coef"]) * highest_exp
                    }
                )"""
        #Setting the objective
        mdlutl.model.objective = mdlutl.model.problem.Objective(Zero, direction="min")
        mdlutl.model.objective.set_linear_coefficients(coef)    
        #Solving the LP
        solution = mdlutl.model.optimize()
        fba_obj = self.save_solution_as_fba(solution,mdlutl,media,params["fba_output_id"],workspace=params["workspace"],fbamodel_ref=params("fbamodel_id"))
        context = {
            "overview": {
                'Model ID':mdlutl.wsid,
                'Media ID':media.id,
                'Target reaction':params["target_reaction"],
                'Gene knockouts':"; ".join(params["feature_ko_list"]),
                'Reaction knockouts':"; ".join(params["reaction_ko_list"]),
                'Exchange limits':exchange_string,
                'Kinetics coefficient':params["kinetics_coef"],
                'Objective fraction':params["objective_fraction"],
                'prFBA':params["prfba"],
                'Max growth':fba_obj.objective_value
            },
            "reactions": [],
            "exchange": [],
            "interaction": []
        }
        env = jinja2.Environment(
            loader=jinja2.FileSystemLoader(self.module_dir+"/data/"),
            autoescape=jinja2.select_autoescape(['html', 'xml']))
        html = env.get_template("FBAReportTemplate.html").render(context)
        os.makedirs(self.working_dir+"/html", exist_ok=True)
        with open(self.working_dir+"/html/index.html", 'w') as f:
            f.write(html)
        #Saving the output
        if params["save_report_to_kbase"]:
            output = self.save_report_to_kbase()
        if bool(params["return_fba_object"]):
            output["fba_obj"] = fba_obj
        return output

    def build_dataframe_report(self,table,model_objs=None):        
        context = {
            "initial_model":table.iloc[0]["Model"]
        }
        env = jinja2.Environment(
            loader=jinja2.FileSystemLoader(self.module_dir+"/data/"),
            autoescape=jinja2.select_autoescape(['html', 'xml']))
        html = env.get_template("ReportTemplate.html").render(context)
        os.makedirs(self.working_dir+"/html", exist_ok=True)
        if model_objs:
            for model in model_objs:
                msmodrep = MSModelReport(model)
                msmodrep.build_report(self.working_dir+"/html/"+model.wsid+"-recon.html")
                msmodrep.build_multitab_report(self.working_dir+"/html/"+model.wsid+"-full.html") 
        print("Output dir:",self.working_dir+"/html/index.html")
        with open(self.working_dir+"/html/index.html", 'w') as f:
            f.write(html)
        #Creating data table file
        json_str = '{"data":'+table.to_json(orient='records')+'}'
        with open(self.working_dir+"/html/data.json", 'w') as f:
            f.write(json_str)