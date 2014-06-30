#!/bin/bash
set -e
set -x

read -p "Enter folder :" id

for folder in $id; do (

echo "Processing $folder"

cd $folder; mkdir -p splitfiles; mkdir -p mcdespot; cd splitfiles

echo "Splitting SPGR and SSFP"
fslsplit ../SPGR_scan.nii.gz spgr -t
fslsplit ../SSFP_scan.nii.gz ssfp -t
fslmaths spgr0000 -subsamp2 target

echo "Registering SPGR images"
for spgr in spgr*; do (
flirt -dof 6 -cost mutualinfo -in $spgr -ref target -out flirt_$spgr.nii.gz
);done

echo "Registering IRSPGR"
flirt -dof 6 -cost mutualinfo -in ../IRSPGR.nii.gz -ref target.nii.gz -out ../mcdespot/irspgr.nii.gz

echo "Registering SSFP to SPGR"
flirt -dof 6 -cost mutualinfo -in ssfp0000.nii.gz -ref target -out flirt_ssfp0000.nii.gz

echo "Registering SSFP images"
for ssfp in ssfp*; do (if [ $ssfp != ssfp0000.nii.gz ] ; then (
flirt -dof 6 -cost mutualinfo -in $ssfp -ref flirt_ssfp0000 -out flirt_$ssfp.nii.gz
); fi ); done

echo "Merging SPGR"
fslmerge -t ../mcdespot/spgr.nii.gz flirt_spgr*.nii.gz

echo "Merging SSFP"
fslmerge -t ../mcdespot/ssfp_both.nii.gz flirt_ssfp*.nii.gz

echo "Creating brain mask"
bet target.nii.gz ../mcdespot/brain -R -m -n

); done
