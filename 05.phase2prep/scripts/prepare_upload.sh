#!/bin/bash

cd ~/sandpit/splitfiles/

tar cvf phase2_list.tar assoclist*
md5sum phase2_list.tar > phase2_list.tar.md5sum
md5sum -c phase2_list.tar.md5sum
