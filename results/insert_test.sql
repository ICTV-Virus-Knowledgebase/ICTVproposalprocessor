--
-- TEST merge_proposals_zip.R output
--
-- SQL to load new MSL, generated directly from proposals.
--

begin transaction

-- TOC insert from OUTPUT
insert into [taxonomy_toc] ([tree_id],[msl_release_num],[comments])  values ( 202200000,38,NULL )

-- QC query
select * from taxonomy_toc
where  msl_release_num =38
select * from taxonomy_node 
where msl_release_num =38

-- TAXONOMY insert from OUTPUT
insert into [taxonomy_node] ([taxnode_id],[parent_id],[tree_id],[msl_release_num],[level_id],[name],[ictv_id],[molecule_id],[abbrev_csv],[genbank_accession_csv],[genbank_refseq_accession_csv],[refseq_accession_csv],[isolate_csv],[notes],[is_hidden],[in_change],[in_target],[in_filename],[in_notes],[out_change],[out_target],[out_filename],[out_notes]) values ('202207117','202200000','202200000','38','120','Duplodnaviria','201907117','1',NULL,NULL,NULL,NULL,NULL,NULL,'1',NULL,NULL,NULL,NULL,NULL,NULL,NULL,NULL) 
insert into [taxonomy_node] ([taxnode_id],[parent_id],[tree_id],[msl_release_num],[level_id],[name],[ictv_id],[molecule_id],[abbrev_csv],[genbank_accession_csv],[genbank_refseq_accession_csv],[refseq_accession_csv],[isolate_csv],[notes],[is_hidden],[in_change],[in_target],[in_filename],[in_notes],[out_change],[out_target],[out_filename],[out_notes]) values ('202207161','202200000','202200000','38','120','Monodnaviria','201907161','2',NULL,NULL,NULL,NULL,NULL,NULL,'1',NULL,NULL,NULL,NULL,NULL,NULL,NULL,NULL) 
insert into [taxonomy_node] ([taxnode_id],[parent_id],[tree_id],[msl_release_num],[level_id],[name],[ictv_id],[molecule_id],[abbrev_csv],[genbank_accession_csv],[genbank_refseq_accession_csv],[refseq_accession_csv],[isolate_csv],[notes],[is_hidden],[in_change],[in_target],[in_filename],[in_notes],[out_change],[out_target],[out_filename],[out_notes]) values ('202207095','202200000','202200000','38','120','Riboviria','201857095',NULL,NULL,NULL,NULL,NULL,NULL,NULL,'1',NULL,NULL,NULL,NULL,NULL,NULL,NULL,NULL) 
insert into [taxonomy_node] ([taxnode_id],[parent_id],[tree_id],[msl_release_num],[level_id],[name],[ictv_id],[molecule_id],[abbrev_csv],[genbank_accession_csv],[genbank_refseq_accession_csv],[refseq_accession_csv],[isolate_csv],[notes],[is_hidden],[in_change],[in_target],[in_filename],[in_notes],[out_change],[out_target],[out_filename],[out_notes]) values ('202209292','202200000','202200000','38','120','Ribozyviria','202009292','5',NULL,NULL,NULL,NULL,NULL,NULL,'1',NULL,NULL,NULL,NULL,NULL,NULL,NULL,NULL) 
insert into [taxonomy_node] ([taxnode_id],[parent_id],[tree_id],[msl_release_num],[level_id],[name],[ictv_id],[molecule_id],[abbrev_csv],[genbank_accession_csv],[genbank_refseq_accession_csv],[refseq_accession_csv],[isolate_csv],[notes],[is_hidden],[in_change],[in_target],[in_filename],[in_notes],[out_change],[out_target],[out_filename],[out_notes]) values ('202208702','202200000','202200000','38','120','Varidnaviria','201908702','1',NULL,NULL,NULL,NULL,NULL,NULL,'1',NULL,NULL,NULL,NULL,NULL,NULL,NULL,NULL) 
insert into [taxonomy_node] ([taxnode_id],[parent_id],[tree_id],[msl_release_num],[level_id],[name],[ictv_id],[molecule_id],[abbrev_csv],[genbank_accession_csv],[genbank_refseq_accession_csv],[refseq_accession_csv],[isolate_csv],[notes],[is_hidden],[in_change],[in_target],[in_filename],[in_notes],[out_change],[out_target],[out_filename],[out_notes]) values ('202207118','202207117','202200000','38','140','Heunggongvirae','201907118',NULL,NULL,NULL,NULL,NULL,NULL,NULL,'1',NULL,NULL,NULL,NULL,NULL,NULL,NULL,NULL) 

--QC query
select * from taxonomy_node 
where msl_release_num =38


rollback transaction
