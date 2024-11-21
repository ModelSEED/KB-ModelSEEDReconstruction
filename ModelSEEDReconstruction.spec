/*
A KBase module: ModelSEEDReconstruction
*/

module ModelSEEDReconstruction {
    typedef structure {
        string report_name;
        string report_ref;
        string workspace;
    } ReportResults;

    /*
        Function builds metabolic models from input genomes
    */
    funcdef build_metabolic_models(mapping<string,UnspecifiedObject> params) returns (ReportResults output) authentication required;
    
    /*
        Function gapfills metabolic models from input genomes
    */
    funcdef gapfill_metabolic_models(mapping<string,UnspecifiedObject> params) returns (ReportResults output) authentication required;

};
