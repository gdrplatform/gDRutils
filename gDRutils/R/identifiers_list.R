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
  drugname = "DrugName",
  drug_moa = "drug_moa",
  # corresponds to the field 'gcsi_drug_name' from gCellGenomics::getDrugs()

  untreated_tag = c("untreated", "vehicle"),
  # flag to identify control treatments

  masked_tag = "masked",
  # flag for masked wells

  well_position = c("WellRow", "WellColumn"),
  concentration = "Concentration",
  template = "Template",
  barcode = "Barcode",
 
  # ids for the 2nd drug 
  drug2 = "Gnumber_2",
  drugname2 = "DrugName_2",
  drug_moa2 = "drug_moa_2",
  concentration2 = "Concentration_2",
  
  # ids for the 3rd drug 
  drug3 = "Gnumber_3",
  drugname3 = "DrugName_3",
  drug_moa3 = "drug_moa_3",
  concentration3 = "Concentration_3"
)
