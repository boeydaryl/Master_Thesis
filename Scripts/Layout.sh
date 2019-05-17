#Bash script to identify Salmon quant folders and save to csv as layout file

#!/bin/bash


for element in */

do

var1=$element
replace=""
var2=${var1///$replace}

echo $var2 >> layout.csv


done