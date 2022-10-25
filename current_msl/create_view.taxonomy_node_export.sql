USE [ICTVonline37]
GO

/****** Object:  View [dbo].[taxonomy_node_export]    Script Date: 10/20/2022 10:29:58 PM ******/
SET ANSI_NULLS ON
GO

SET QUOTED_IDENTIFIER ON
GO

CREATE view [dbo].[taxonomy_node_export] as
--
-- Export used by validate_proposals.Rmd load script
-- Should be saved as current_msl/taxonony_node_export.txt
-- As that will have no column headers
--
-- export: select * from taxonomy_node_export 
select 
	top 10000000 -- add top so we can use 'order by' in a view
     tn.[taxnode_id]
    ,tn.[parent_id]
    ,tn.[tree_id]
    ,tn.[msl_release_num]
    ,tn.[level_id]
    ,tn.[name]
    ,tn.[ictv_id]
    ,tn.[molecule_id]
    ,tn.[abbrev_csv]
    ,tn.[genbank_accession_csv]
    ,tn.[genbank_refseq_accession_csv]
    ,tn.[refseq_accession_csv]
    ,tn.[isolate_csv]
    ,tn.[notes]
    ,tn.[is_ref]
    ,tn.[is_official]
    ,tn.[is_hidden]
    ,tn.[is_deleted]
    ,tn.[is_deleted_next_year]
    ,tn.[is_typo]
    ,tn.[is_renamed_next_year]
    ,tn.[is_obsolete]
    ,tn.[in_change]
    ,tn.[in_target]
    ,tn.[in_filename]
    ,tn.[in_notes]
    ,tn.[out_change]
    ,tn.[out_target]
    ,tn.[out_filename]
    ,tn.[out_notes]
	-- Trigger-maintained and computed columns
	,tn.[lineage]
	,tn.[cleaned_name]
	-- LUT/join values
	,[rank]   = isnull([rank].name, '')
	,[molecule] = isnull([mol].abbrev,'')
from taxonomy_node tn
left outer join taxonomy_level [rank] on [rank].id=tn.level_id
left outer join taxonomy_molecule mol on mol.id = tn.molecule_id
where tn.msl_release_num is not null
and tn.is_deleted = 0 
and (tn.level_id = 100 or tn.is_hidden = 0) 
and tn.is_obsolete=0
order by tn.msl_release_num, tn.left_idx
GO

