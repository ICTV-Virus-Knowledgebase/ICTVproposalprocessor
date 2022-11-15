--
-- TEST merge_proposals_zip.R output
--
-- SQL to load new MSL, generated directly from proposals.
--
-- https://github.com/ICTV-Virus-Knowledgebase/MSL_merge/blob/check/results/insert_test.sql


begin transaction

-- QC query
select [taxonomy_toc.before]='before', * from taxonomy_toc
where  msl_release_num =38
select [taxonomy_node.before]='before',* from taxonomy_node 
where msl_release_num =38

-- TAXONOMY insert from OUTPUT




--QC query
select [taxonomy_toc.after]='after', * from taxonomy_toc
where  msl_release_num =38

select  [taxonomy_node.after]='after',* from taxonomy_node 
where msl_release_num =38
order by left_idx



rollback transaction
