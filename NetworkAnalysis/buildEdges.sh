
#!/usr/local/bin/bash

## $1 is the pdb file (X.pdb)

out=${1%%.pdb}

arr=($(echo $out | tr "/" "\n"))
dir=${arr[0]}
file=${arr[1]}

echo "Writing out polypeptide file..."
python pdb2polypeptide.py $1
awk '{print $1"\t"NR; print$2"\t"NR+1}' "$out"_polypeptidePairs | uniq > "$out"_AAorder

echo "Adding hydrogens..."
#Ensure that this line reflects the path of your Phenix installation.
source ../../../../../../../../Applications/phenix-1.14-3260/phenix_env.sh 
phenix.reduce -NOFLIP -Quiet $out".pdb" > $out"_H.pdb"
mv $out"_H.pdb" $out".pdb"

echo "Protonating waters..."
if grep -q "HOH" $out".pdb"
then
sed s/"AHOH"/" HOH"/g $out".pdb" | sed /"BHOH"/d > $out"_tmp.pdb"
phenix.ready_set add_h_to_water=True output_file_name=$out"_H" $out"_tmp.pdb"
rm $out"_tmp.pdb"
mv $out"_H.pdb" $out".pdb"
fi

echo "Making .phi file"
../Stride/stride $out".pdb" > $out".phi"

echo "Running full bond calculation script"
python pdb2edgeSC.py $out".pdb"

echo "Finding water bonds"
python makeWaterBonds.py $out"_hb"
awk '{split($1,x,"-");split($2,y,"-");if(x[1]x[2]!=y[1]y[2]){print}}' $out"_waterBonds" > $out"tmp"
mv $out"tmp" $out"_waterBonds"
    
echo "Making "$out"_net"
awk '{print $1"\t"$2"\tMCMC\t10\tPP\tPP1\tPP2"}' $out"_polypeptidePairs" > $out"_net"
awk '{split($1,x,"-");if(s[x[1]"\t"x[2]]==""&&s[x[2]"\t"x[1]]==""){s[x[1]"\t"x[2]]=$2;next;}if(s[x[1]"\t"x[2]]!=""&&s[x[2]"\t"x[1]]==""){s[x[1]"\t"x[2]]=s[x[1]"\t"x[2]]+$2;next;}if(s[x[1]"\t"x[2]]==""&&s[x[2]"\t"x[1]]!=""){s[x[2]"\t"x[1]]=s[x[2]"\t"x[1]]+$2}}END{for(i in s){if(s[i]!=""&&s[i]>0){print i"\t"s[i]}}}' $out"_vdw" > $out"tmp"
awk '{print $1"-"$2;print $2"-"$1;}' $out"tmp" > $out"tmp2"
grep -f $out"tmp2" $out"_vdw" > $out"_vdw_noRepulse"
awk '{split($1,x,"-");if((x[1]~/PRO/&&x[5]=="CD")||(x[2]~/PRO/&&x[6]=="CD")){print x[1]"\t"x[2]"\t"x[3]x[4]"\t"$2+0"\tVDW\t"$3"\t"$4;}else{print x[1]"\t"x[2]"\t"x[3]x[4]"\t"$2"\tVDW\t"$3"\t"$4;}}' $out"_vdw_noRepulse" >> $out"_net"
awk '{print $1"\t"$2"\tSCSC\t"$3"\tPIPI\t"$4"\t"$5}' $out"_pipi2" >> $out"_net"
awk '{print $1"\t"$2"\tSCSC\t"$3"\tPICAT\t"$4"\t"$5}' $out"_pication2" >> $out"_net"
awk '{print $1"\t"$2"\tSCSC\t"$3"\tSS\t"$4"\t"$5}' $out"_disulf" >> $out"_net"
awk '{split($1,x,"-");split($2,y,"-");print x[1]x[2]"\t"y[1]y[2]"\tSCSC\t"$3"\tSB\t"$4"\t"$5}' $out"_saltBridges_Barlow" >> $out"_net"
awk '{split($1,x,"-");split($2,y,"-");print x[1]x[2]"\t"y[1]y[2]"\t"$3"\t"$10"\tHB\t"$11"\t"$12}' $out"_hb" | sed /HOH/d >> $out"_net"
awk '{split($1,x,"-");split($2,y,"-");print x[1]x[2]"\t"y[1]y[2]"\t"$3"\t"$4"\t"$5"\t"$6"\t"$7}' $out"_waterBonds" >> $out"_net"
awk '{split($1,x,"-");split($2,y,"-");print x[1]x[2]"\t"y[1]y[2]"\t"$3"\t"3"\tMETAL\t"x[3]"\t"y[3]}' $out"_metal" >> $out"_net"
awk '{split($1,x,"-");split($2,y,"-");print x[1]x[2]"\t"y[1]y[2]"\tMC"$3"\t10\tDNA\tNT\t"x[3]}' $out"_DNA" >> $out"_net"


awk '{split($1,x,"-");if(s[x[1]"\t"x[2]]==""&&s[x[2]"\t"x[1]]==""){s[x[1]"\t"x[2]]=$2;next;}if(s[x[1]"\t"x[2]]!=""&&s[x[2]"\t"x[1]]==""){s[x[1]"\t"x[2]]=s[x[1]"\t"x[2]]+$2;next;}if(s[x[1]"\t"x[2]]==""&&s[x[2]"\t"x[1]]!=""){s[x[2]"\t"x[1]]=s[x[2]"\t"x[1]]+$2}}END{for(i in s){if(s[i]!=""&&s[i]>0){print i"\t"s[i]}}}' $out"_vdw2" > $out"tmp"
awk '{print $1"-"$2;print $2"-"$1;}' $out"tmp" > $out"tmp2"
grep -f $out"tmp2" $out"_vdw2" > $out"_vdw2_noRepulse"
awk '{split($1,x,"-");if((x[1]~/PRO/&&x[5]=="CD")||(x[2]~/PRO/&&x[6]=="CD")){print x[1]"\t"x[2]"\t"x[3]x[4]"\t"$2+0"\tVDW2\t"$3"\t"$4;}else{print x[1]"\t"x[2]"\t"x[3]x[4]"\t"$2"\tVDW2\t"$3"\t"$4;}}' $out"_vdw2_noRepulse" >> $out"_net"

echo "Removing any negative value edges"
awk '$4>0{print}' $out"_net" > $out"tmp"
mv $out"tmp" $out"_net"

echo "Creating QC file"
wc -l $out"_"* > $out".QC"

echo "Removing water from centroid"
sed /HOH/d $out"_centroidNetSC" > tmp
mv tmp $out"_centroidNetSC"


cd $dir
mkdir Centroid
cp $file* Centroid
cd Centroid
mv $file"_centroidNetSC" $file"_net"
cd ../..

##Run network analysis
R --vanilla --args $out < modularity_analysis_energetics.R
cd $dir
cd Centroid
R --vanilla --args $file < ../../modularity_analysis_centroid.R
