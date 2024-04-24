#!/usr/bin/env bash
#
# list proposal names
#

find . -type f | egrep -v "(.DS_Store|file_list|downloads)" > file_list.txt
find . -type f -name "20*" -exec basename "{}" .x \; | egrep -v "(.DS_Store|file_list|downloads)" > file_names.txt

