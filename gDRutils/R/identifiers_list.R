## The following are the default values that will be used for critical metadata fields
## if no overriding is performed via `set_identifier()`.

IDENTIFIERS_LIST <- list(
  duration = "Duration",

  cellline = "clid",
  cellline_name = "CellLineName",
  cellline_tissue = "Tissue",
  cellline_ref_div_time = "ReferenceDivisionTime",

  drug = "Gnumber",
  drugname = "DrugName",
  # corresponds to the field 'gcsi_drug_name' from gCellGenomics::getDrugs()

  untreated_tag = c("untreated", "vehicle"),
  # flag to identify control treatments

  masked_tag = 'masked',
  # flag for masked wells

  WellPosition = c("WellRow", "WellColumn"), 
  # TODO: Delete once all references to WellPosition are replaced with well_position
  well_position = c("WellRow", "WellColumn")
)
