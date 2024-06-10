#!/usr/bin/env bash
#
# compare pre- and post- name cleanup lists
#
# from finalize_proposals.R
#
DIR1=Pending_Proposals
DIR2=proposalsFinal

# create clean file lists
echo "build file_list: $DIR1"
(cd $DIR1; ../file_list.sh)
echo "build file_list: $DIR2"
(cd $DIR2; ../file_list.sh)


# compare
echo "Compare $DIR1 vs $DIR2"
diff \
    <(cat $DIR1/file_names.txt| sed -e 's/\.A\././;s/\.v[0-9]\././'|sort) \
    <(cat $DIR2/file_names.txt| grep -v download |sort)
