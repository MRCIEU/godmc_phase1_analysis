#!/bin/bash
ls

status=$?
if [ ! "$status" -eq "0" ]
then
	echo $status
fi
echo $status

