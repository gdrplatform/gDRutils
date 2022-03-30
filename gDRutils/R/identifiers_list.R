## The following are the default values that will be used for critical metadata fields
## if no overriding is performed via `set_identifier()`.

IDENTIFIERS_LIST <- list(
  duration = "Duration",

  cellline = "clid",
  cellline_name = "CellLineName",
  cellline_tissue = "Tissue",
  cellline_ref_div_time = "ReferenceDivisionTime",
  cellline_parental_identifier = "parental_identifier",
  cellline_subtype = "subtype",

  drug = "Gnumber",
  drug_name = "DrugName",
  drug_moa = "drug_moa",
  # corresponds to the field 'gcsi_drug_name' from gCellGenomics::getDrugs()

  untreated_tag = c("untreated", "vehicle"),
  # flag to identify control treatments

  masked_tag = "masked",
  # flag for masked wells

  well_position = c("WellRow", "WellColumn"),
  concentration = "Concentration",
  template = "Template",
  barcode = c("Barcode", "Plate"),
 
  # ids for the 2nd drug 
  drug2 = "Gnumber_2",
  drug_name2 = "DrugName_2",
  drug_moa2 = "drug_moa_2",
  concentration2 = "Concentration_2",
  
  # ids for the 3rd drug 
  drug3 = "Gnumber_3",
  drug_name3 = "DrugName_3",
  drug_moa3 = "drug_moa_3",
  concentration3 = "Concentration_3",
  
  # data source
  data_source = "data_source"
)

REQ_COL_IDENTIFIERS <- c(
  "duration",
  "cellline",
  "cellline_name",
  "drug",
  "drug_name",
  "masked_tag",
  "concentration"
)

EXPECT_ONE_IDENTIFIERS <- c(
  "duration",
  "cellline",
  "cellline_name",
  "cellline_tissue",
  "cellline_ref_div_time",
  "cellline_parental_identifier",
  "cellline_subtype",
  "drug",
  "drug_name",
  "drug_moa",
  "masked_tag",
  "concentration",
  "template",
  "barcode",
  "drug2",
  "drug_name2",
  "drug_moa2",
  "concentration2",
  "drug3",
  "drug_name3",
  "drug_moa3",
  "concentration3",
  "data_source"
)
