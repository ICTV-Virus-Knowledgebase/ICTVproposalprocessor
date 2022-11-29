--
-- TEST merge_proposals_zip.R output
--
-- SQL to load new MSL, generated directly from proposals.
--
-- https://github.com/ICTV-Virus-Knowledgebase/MSL_merge/blob/check/results/insert_test.sql

begin transaction
-- rollback transaction

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

-- report on changes in taxonomy_node
select msl_release_num, level_id, change='in_change', in_ct=count(in_change), ct=count(*)
  from taxonomy_node
  where msl_release_num = 38
  group by msl_release_num, level_id
union all
select msl_release_num, level_id, change='out_change', col_ct=count(out_change), ct=count(*)
  from taxonomy_node
  where msl_release_num = 37
  group by msl_release_num, level_id
order by msl_release_num, level_id


exec rebuild_delta_nodes NULL

-- report on deltas
select prev_tags,ct=count(*) from taxonomy_node_dx where msl_release_num = 38  group by prev_tags


rollback transaction
