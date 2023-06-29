#!/bin/bash

# #####################################################################################
# # Script functions and settings

shopt -s nullglob

# Create our "press_enter function"
function press_enter
{
    echo ""
    echo -n "Press Enter to continue"
    read
    clear
}

# Create a simple method for using comment blocks
[ -z $BASH ] || shopt -s expand_aliases
alias BCOMM="if [ ]; then"
alias ECOMM="fi"

# #####################################################################################

# if [ "$1" == "-h" ]; then
#   echo "------------------------------------------------------------------------------------------------------------------------"
#   echo "Pipeline for image analysis for CALIPR framework study"
#   echo ""
#   echo " Run inside of CALIPR data directory"
#   echo ""
#   echo " Written by Adam Dvorak (2023)"
#   echo ""
#   echo "------------------------------------------------------------------------------------------------------------------------"
#   exit 0
# fi

###############  SETUP PATHS  ###############
# # TITAN
# dataPath='/home/bizon/Research/MWI_Development/Data'
# inputPath='/home/bizon/Research/MWI_Development/CALIPR/CALIPR_Repro'
# oasisPath='/home/bizon/Research/Templates/OASIS'
# qcPath=${inputPath}/QualityControl
# cores=30
# export ITK_GLOBAL_DEFAULT_NUMBER_OF_THREADS=30


# MBP
dataPath='/Users/adamdvorak/Research/MWI_Development/Data'
inputPath='/Users/adamdvorak/Research/MWI_Development/CALIPR/MWI_Repro'
oasisPath='/Users/adamdvorak/Research/Scripts_Tools/Atlases/OASIS/'
qcPath=${inputPath}/QualityControl
grase48Path='/Users/adamdvorak/Research/Archived/Atlas/Template/'
cores=8
export ITK_GLOBAL_DEFAULT_NUMBER_OF_THREADS=8
###############  SETUP PATHS  ###############


###############  SETUP SUBJECTS  ###############
subjects=' calipr_repro_c01 calipr_repro_c02 calipr_repro_c03 calipr_repro_c04 calipr_repro_c05 '
# subjects=' lobstr_bisque_pilot '
###############  SETUP SUBJECTS  ###############

# always start with cd into inputPath
cd ${inputPath}



#######################################################################################################################################################################################
######################################################################################## Cord ########################################################################################
#######################################################################################################################################################################################


# BCOMM
##################################################################################### Cord Prep, MWI<->mFFE Registration #####################################################################################
for subject in ${subjects}
do
  cd ${inputPath}/${subject}/cord
  printf " \n Beginning ${subject} Cord Prep \n "

  ##### mFFE
  sct_propseg \
    -i ${subject}_cord_mffe.nii.gz \
    -c t2s -CSF -qc ${qcPath}

  sct_deepseg_gm \
    -i ${subject}_cord_mffe.nii.gz \
    -m large -thr 0.999 -qc ${qcPath}

  mv ${subject}_cord_mffe_gmseg.nii.gz ${subject}_cord_mffe_segGM.nii.gz
  mv ${subject}_cord_mffe_CSF_seg.nii.gz ${subject}_cord_mffe_segCSF.nii.gz

  # make combined cord+csf seg
  fslmaths ${subject}_cord_mffe_seg.nii.gz -add ${subject}_cord_mffe_segCSF.nii.gz -bin -fillh ${subject}_cord_mffe_segwCSF.nii.gz
  # make wm seg
  fslmaths ${subject}_cord_mffe_seg.nii.gz -sub ${subject}_cord_mffe_segGM.nii.gz -thr 1 -bin ${subject}_cord_mffe_segWM.nii.gz
  # mask
  fslmaths ${subject}_cord_mffe.nii.gz -mas ${subject}_cord_mffe_segwCSF.nii.gz ${subject}_cord_mffe_cordwCSF.nii.gz
  fslmaths ${subject}_cord_mffe.nii.gz -mas ${subject}_cord_mffe_seg.nii.gz ${subject}_cord_mffe_cord.nii.gz

  ##### MWI
  for mwi in $(ls  *cord_calipr_1.nii* *cord_calipr_2.nii* )
  do
    printf " \n Beginning ${subject} ${mwi%.*} Preparation \n "

    fslroi ${mwi%.*}.nii ${mwi%.*}_E1.nii.gz 0 1
    fslroi ${mwi%.*}.nii ${mwi%.*}_E48.nii.gz 47 1

    sct_propseg -i ${mwi%.*}_E48.nii.gz -c t2 -CSF -qc ${qcPath}

    mv ${mwi%.*}_E48_CSF_seg.nii.gz ${mwi%.*}_E48_segCSF.nii.gz

    # make combined cord+csf seg
    fslmaths ${mwi%.*}_E48_seg.nii.gz -add ${mwi%.*}_E48_segCSF.nii.gz -bin -fillh ${mwi%.*}_E48_segwCSF.nii.gz
    # mask
    fslmaths ${mwi%.*}_E48.nii.gz -mas ${mwi%.*}_E48_segwCSF.nii.gz ${mwi%.*}_E48_cordwCSF.nii.gz
    fslmaths ${mwi%.*}_E48.nii.gz -mas ${mwi%.*}_E48_seg.nii.gz ${mwi%.*}_E48_cord.nii.gz

    # Register
      antsRegistrationSyN.sh \
      -d 3 \
      -f ${subject}_cord_mffe_cordwCSF.nii.gz \
      -m ${mwi%.*}_E48_cordwCSF.nii.gz \
      -t r \
      -z 1 \
      -j 0 \
      -p d \
      -n ${cores} \
      -o ${mwi%.*}_E48_cordwCSF
  done
done
##################################################################################### Cord Prep, MWI<->mFFE Registration #####################################################################################
# ECOMM


# BCOMM
##################################################################################### Warp MWI to mFFE #####################################################################################
for subject in ${subjects}
do
  cd ${inputPath}/${subject}/cord
  printf " \n Beginning ${subject} Warp MWI to mFFE \n "

  ##### MWI
  for mwi in $(ls  *cord_calipr_1.nii* *cord_calipr_2.nii* )
  do

    ########## Warp full volume, and create N4 correcting version
    # Redo N4 correction, saving bias field
    N4BiasFieldCorrection \
      -d 3 \
      -i ${mwi%.*}_E1.nii.gz \
      -o [ ${mwi%.*}_E1_N4.nii.gz, ${mwi%.*}_E1_N4biasfield.nii.gz ]

    # make N4 corrected version
    fslmaths \
      ${mwi} \
      -div ${mwi%.*}_E1_N4biasfield.nii.gz \
      ${mwi%.*}_N4.nii.gz
    # warp
    antsApplyTransforms \
      -d 3 \
      -e 3 \
      -i ${mwi%.*}_N4.nii.gz \
      -r ${subject}_cord_mffe.nii.gz \
      -t ${mwi%.*}_E48_cordwCSF0GenericAffine.mat \
      -o ${mwi%.*}_N4_mffeWarped.nii.gz


    ##### images/masks
    towarp=' E1 E48 E48_seg E48_segwCSF '
    for image in ${towarp}
    do

      # ls ${mwi%.*}_${image}.nii*

      antsApplyTransforms \
        -d 3 \
        -i ${mwi%.*}_${image}.nii* \
        -r ${subject}_cord_mffe.nii.gz \
        -t ${mwi%.*}_E48_cordwCSF0GenericAffine.mat \
        -o ${mwi%.*}_${image}_mffeWarped.nii.gz

    done

    # Re-threshold warped segmentations
    fslmaths ${mwi%.*}_E48_seg_mffeWarped.nii.gz -thr 0.99 -bin ${mwi%.*}_E48_seg_mffeWarped.nii.gz
    fslmaths ${mwi%.*}_E48_segwCSF_mffeWarped.nii.gz -thr 0.99 -bin ${mwi%.*}_E48_segwCSF_mffeWarped.nii.gz

    ##### quantitative maps
    for iter in {0..2}
    do
      towarp=' MWF IET2 ALPHA '
      for image in ${towarp}
      do

        # ls ${mwi%.*}_DK${iter}_${image}.nii*

        antsApplyTransforms \
          -d 3 \
          -i ${mwi%.*}_DK${iter}_${image}.nii* \
          -r ${subject}_cord_mffe.nii.gz \
          -t ${mwi%.*}_E48_cordwCSF0GenericAffine.mat \
          -o ${mwi%.*}_DK${iter}_${image}_mffeWarped.nii.gz

      done
    done

  done
done
##################################################################################### Warp MWI to mFFE #####################################################################################
# ECOMM



# BCOMM
##################################################################################### mFFE<->Template Reg #####################################################################################
for subject in ${subjects}
do
  cd ${inputPath}/${subject}/cord
  printf " \n Beginning ${subject} Cord mFFE<->Template Reg \n "

  # centred at level of C3/C4 disc
  sct_label_utils \
    -i ${subject}_cord_mffe_seg.nii.gz \
    -create-seg-mid 4 \
    -o ${subject}_cord_mffe_seg_labeled_disc.nii.gz

  sct_register_to_template \
    -i ${subject}_cord_mffe.nii.gz \
    -s ${subject}_cord_mffe_segGM.nii.gz \
    -s-template-id 5 \
    -ldisc ${subject}_cord_mffe_seg_labeled_disc.nii.gz \
    -c t2s \
    -ref subject \
    -ofolder mffeTemplateReg

  sct_warp_template \
    -d ${subject}_cord_mffe.nii.gz \
    -w mffeTemplateReg/warp_template2anat.nii.gz \
    -ofolder mffeTemplateReg

done
##################################################################################### mFFE<->Template Reg #####################################################################################
# ECOMM



# BCOMM
############################################################################## Cord mask and ROI Prep ##############################################################################
for subject in ${subjects}
do
  
  # Change into the subject folder
  cd ${inputPath}/${subject}/cord

  ############################################################ Create intitial ROIs for extraction
  fslmerge -t mffeTemplateReg/atlas/PAM50_atlas_LF.nii.gz $( for i in $(seq -f "%02g" 4 13); do printf "mffeTemplateReg/atlas/PAM50_atlas_${i}.nii.gz "; done )
  fslmaths mffeTemplateReg/atlas/PAM50_atlas_LF.nii.gz -Tmax -thr 0.5 -bin mffeTemplateReg/atlas/PAM50_atlas_LF.nii.gz

  fslmerge -t mffeTemplateReg/atlas/PAM50_atlas_VF.nii.gz $( for i in $(seq -f "%02g" 14 29); do printf "mffeTemplateReg/atlas/PAM50_atlas_${i}.nii.gz "; done )
  fslmaths mffeTemplateReg/atlas/PAM50_atlas_VF.nii.gz -Tmax -thr 0.5 -bin mffeTemplateReg/atlas/PAM50_atlas_VF.nii.gz

  fslmerge -t mffeTemplateReg/atlas/PAM50_atlas_DC.nii.gz $( for i in $(seq -f "%02g" 0 3); do printf "mffeTemplateReg/atlas/PAM50_atlas_${i}.nii.gz "; done )
  fslmaths mffeTemplateReg/atlas/PAM50_atlas_DC.nii.gz -Tmax -thr 0.5 -bin mffeTemplateReg/atlas/PAM50_atlas_DC.nii.gz

  fslmerge -t mffeTemplateReg/atlas/PAM50_atlas_LCST.nii.gz $( for i in $(seq -f "%02g" 4 5); do printf "mffeTemplateReg/atlas/PAM50_atlas_${i}.nii.gz "; done )
  fslmaths mffeTemplateReg/atlas/PAM50_atlas_LCST.nii.gz -Tmax -thr 0.5 -bin mffeTemplateReg/atlas/PAM50_atlas_LCST.nii.gz

  # Mask and prep final ROI versions:
  mkdir ROIs

  fslmaths mffeTemplateReg/template/PAM50_cord.nii.gz -thr 0.5 -bin -mas ${subject}_cord_mffe_seg.nii.gz ROIs/ROI_WC.nii.gz
  fslmaths mffeTemplateReg/template/PAM50_gm.nii.gz -thr 0.5 -bin -mas ${subject}_cord_mffe_segGM.nii.gz -mas ROIs/ROI_WC.nii.gz ROIs/ROI_GM.nii.gz
  fslmaths mffeTemplateReg/template/PAM50_wm.nii.gz -thr 0.5 -bin -mas ${subject}_cord_mffe_segWM.nii.gz -mas ROIs/ROI_WC.nii.gz ROIs/ROI_WM.nii.gz


  rois=' LF VF DC LCST '
  for roi in ${rois}
  do
    fslmaths mffeTemplateReg/atlas/PAM50_atlas_${roi}.nii.gz -mas ROIs/ROI_WM.nii.gz ROIs/ROI_${roi}.nii.gz
  done


  ################### Create MWI Overlap Masks
  echo "Creating MWI overlap mask from:"
  ls *cord_calipr*_seg_mffeWarped.nii.gz

  fsladd \
    ${subject}_cord_mffe_MaskMWIoverlap.nii.gz -m \
    *cord_calipr*_seg_mffeWarped.nii.gz
  
  fslmaths \
    ${subject}_cord_mffe_MaskMWIoverlap.nii.gz -thr 1 \
    ${subject}_cord_mffe_MaskMWIoverlap.nii.gz


  # finally, mask all ROIs to MWIoverlap
  rois=' WC WM GM LF VF DC LCST '
  for roi in ${rois}
  do
    fslmaths ROIs/ROI_${roi}.nii.gz -mas ${subject}_cord_mffe_MaskMWIoverlap.nii.gz ROIs/ROI_${roi}_Masked.nii.gz
  done

  printf " \n ${subject} ROI prep complete \n "

done
############################################################################## Cord mask and ROI Prep ##############################################################################
# ECOMM



# BCOMM
######################################################################## Create MWI Difference Maps ########################################################################
for subject in ${subjects}
do
  # Change into the subject folder
  cd ${inputPath}/${subject}/cord

  ##### quantitative maps
  for iter in {0..2}
  do
    towarp=' MWF IET2 '
    for image in ${towarp}
    do

      # mask with sharper 3DT1 brain mask
      fslmaths \
        ${subject}_cord_calipr_1_DK${iter}_${image}_mffeWarped.nii.gz \
        -sub ${subject}_cord_calipr_2_DK${iter}_${image}_mffeWarped.nii.gz \
        ${subject}_cord_calipr_1sub2_DK${iter}_${image}_mffeWarped.nii.gz
      # absolute value version
      fslmaths \
        ${subject}_cord_calipr_1sub2_DK${iter}_${image}_mffeWarped.nii.gz \
        -abs ${subject}_cord_calipr_1sub2abs_DK${iter}_${image}_mffeWarped.nii.gz

      # cord-masked, MWI-overlap-masked versions of each
      fslmaths \
        ${subject}_cord_calipr_1sub2_DK${iter}_${image}_mffeWarped.nii.gz \
        -mas ROIs/ROI_WC_Masked.nii.gz \
        ${subject}_cord_calipr_1sub2_DK${iter}_${image}_mffeWarpedMasked.nii.gz
      fslmaths \
        ${subject}_cord_calipr_1sub2abs_DK${iter}_${image}_mffeWarped.nii.gz \
        -mas ROIs/ROI_WC_Masked.nii.gz \
        ${subject}_cord_calipr_1sub2abs_DK${iter}_${image}_mffeWarpedMasked.nii.gz

    done
  done
done
######################################################################## Create MWI Difference Maps ########################################################################
# ECOMM


# BCOMM
######################################################################## Create Masked Versions of Warped Maps ########################################################################
for subject in ${subjects}
do
  cd ${inputPath}/${subject}/cord
  printf " \n ${subject} creating masked versions of warped maps (in mffe space) \n "

  ##### MWI
  for mwi in $(ls  *cord_calipr_1.nii* *cord_calipr_2.nii* )
  do

    ##### images/masks
    towarp=' E1 E48 E48_seg E48_segwCSF '
    for image in ${towarp}
    do

      fslmaths \
        ${mwi%.*}_${image}_mffeWarped.nii.gz \
        -mas ROIs/ROI_WC_Masked.nii.gz \
        ${mwi%.*}_${image}_mffeWarpedMasked.nii.gz

    done

    ##### quantitative maps
    for iter in {0..2}
    do
      towarp=' MWF IET2 ALPHA '
      for image in ${towarp}
      do

        fslmaths \
          ${mwi%.*}_DK${iter}_${image}_mffeWarped.nii.gz \
          -mas ROIs/ROI_WC_Masked.nii.gz \
          ${mwi%.*}_DK${iter}_${image}_mffeWarpedMasked.nii.gz

      done
    done
  done
done
######################################################################## Create Masked Versions of Warped Maps ########################################################################
# ECOMM



############################################################################## EXTRACT ROI Values ##############################################################################


# BCOMM
##################################################################################### Extract for entire ROI
# rois=' WMandGM WM GM JHU CC Genu Splenium PIC Caudate Putamen Thalamus Cortical_GM Frontal_WM Occipital_WM Parietal_WM Temporal_WM '

rois=' WC WM GM LF VF DC LCST '
# rois=' WC WM GM '


rm ${inputPath}/ROIresults.txt
# Create Header
printf "Subject  Data  Metric  ROI  Mask  ROIMask  Mean  Median  Stddev  5perc  95perc  Volume  Volume_mm3  DKiter  \n" >> ${inputPath}/ROIresults.txt

for subject in ${subjects}
do
  for iter in {0..2}
  do
    # Change into the subject folder
    cd ${inputPath}/${subject}/cord

    # for each MWI dataset
    for mwi in $(ls  *cord_calipr_1.nii* *cord_calipr_2.nii* )
    do
      imgs=' MWF IET2 '
      # imgs=' MWF '
      for img in ${imgs}
      do
        
        printf " \n Extracting ${mwi%.*}   DK${iter}   ${img} \n \n "
        ###############################   WMandGM
        for roi in ${rois}
        do
          # loop over roi
          printf "${subject}  ${mwi%.*}  ${img}  ${roi}  Cord  MWIoverlap  " >> ${inputPath}/ROIresults.txt
          # mean
          fslstats ${mwi%.*}_DK${iter}_${img}_mffeWarpedMasked.nii.gz -k ROIs/ROI_${roi}_Masked.nii.gz -m | tr '\n' ' ' >> ${inputPath}/ROIresults.txt
          # median
          fslstats ${mwi%.*}_DK${iter}_${img}_mffeWarpedMasked.nii.gz -k ROIs/ROI_${roi}_Masked.nii.gz -p 50 | tr '\n' ' ' >> ${inputPath}/ROIresults.txt
          # std
          fslstats ${mwi%.*}_DK${iter}_${img}_mffeWarpedMasked.nii.gz -k ROIs/ROI_${roi}_Masked.nii.gz -s | tr '\n' ' ' >> ${inputPath}/ROIresults.txt
          # 5th percential
          fslstats ${mwi%.*}_DK${iter}_${img}_mffeWarpedMasked.nii.gz -k ROIs/ROI_${roi}_Masked.nii.gz -p 5 | tr '\n' ' ' >> ${inputPath}/ROIresults.txt
          # 95th percential
          fslstats ${mwi%.*}_DK${iter}_${img}_mffeWarpedMasked.nii.gz -k ROIs/ROI_${roi}_Masked.nii.gz -p 95 | tr '\n' ' ' >> ${inputPath}/ROIresults.txt
          # size of ROI in this image
          fslstats ROIs/ROI_${roi}_Masked.nii.gz -V | tr '\n' ' ' >> ${inputPath}/ROIresults.txt
          # iterations of DK analysis
          printf "${iter}  \n" >> ${inputPath}/ROIresults.txt
        done
      done
      imgs=' ALPHA '
      for img in ${imgs}
      do
        printf " \n Extracting ${mwi%.*}   DK${iter}   ${img} \n \n "
        ###############################   WMandGM
        for roi in WC
        do
          # loop over roi
          printf "${subject}  ${mwi%.*}  ${img}  ${roi}  Cord  MWIoverlap  " >> ${inputPath}/ROIresults.txt
          # mean
          fslstats ${mwi%.*}_DK${iter}_${img}_mffeWarpedMasked.nii.gz -k ROIs/ROI_${roi}_Masked.nii.gz -m | tr '\n' ' ' >> ${inputPath}/ROIresults.txt
          # median
          fslstats ${mwi%.*}_DK${iter}_${img}_mffeWarpedMasked.nii.gz -k ROIs/ROI_${roi}_Masked.nii.gz -p 50 | tr '\n' ' ' >> ${inputPath}/ROIresults.txt
          # std
          fslstats ${mwi%.*}_DK${iter}_${img}_mffeWarpedMasked.nii.gz -k ROIs/ROI_${roi}_Masked.nii.gz -s | tr '\n' ' ' >> ${inputPath}/ROIresults.txt
          # 5th percential
          fslstats ${mwi%.*}_DK${iter}_${img}_mffeWarpedMasked.nii.gz -k ROIs/ROI_${roi}_Masked.nii.gz -p 5 | tr '\n' ' ' >> ${inputPath}/ROIresults.txt
          # 95th percential
          fslstats ${mwi%.*}_DK${iter}_${img}_mffeWarpedMasked.nii.gz -k ROIs/ROI_${roi}_Masked.nii.gz -p 95 | tr '\n' ' ' >> ${inputPath}/ROIresults.txt
          # size of ROI in this image
          fslstats ROIs/ROI_${roi}_Masked.nii.gz -V | tr '\n' ' ' >> ${inputPath}/ROIresults.txt
          # iterations of DK analysis
          printf "${iter}  \n" >> ${inputPath}/ROIresults.txt
        done
      done
    done
  done
done

#####################################################################################
# ECOMM






#######################################################################################################################################################################################
######################################################################################## BRAIN ########################################################################################
#######################################################################################################################################################################################

# BCOMM
##################################################################################### STRUCTURE PREP #####################################################################################

# rename calipr_cord -> cord_calipr, same for brain
for file in *_calipr_brain_* *_calipr_cord_*
do
  mv "$file" "${file/_calipr_brain_/_brain_calipr_}"
  mv "$file" "${file/_calipr_cord_/_cord_calipr_}"
done

# convert any anatomical par/recs in folder
for file in *_brain_t1w.par *_cord_mffe.par
do
  dcm2niix -b n -s y -m y -p y -z y -f "${file/.par/}" ${file}
done

# clean up all the MWI echoes for now
rm *_e*.nii.gz*

# SORT
for subject in ${subjects}
do
  mkdir ${inputPath}/${subject} ${inputPath}/${subject}/brain ${inputPath}/${subject}/cord

  mv ${inputPath}/${subject}*brain* ${subject}/brain/
  mv ${inputPath}/${subject}*cord* ${subject}/cord/
done

##################################################################################### STRUCTURE PREP #####################################################################################
# ECOMM


# BCOMM
##################################################################################### BRAIN 3DT1 #####################################################################################

cd ${inputPath}
mkdir ${qcPath}

for subject in ${subjects}
do
  cd ${inputPath}/${subject}/brain

  # start the timer
  timer_start="$(date +"Date : %d/%m/%Y Time : %H.%M.%S")"
  printf " \n Beginning ${subject} t1w Preparation at: \n ${timer_start} \n "

  # N4 correction
  N4BiasFieldCorrection \
    -d 3 \
    -i ${subject}_brain_t1w.nii.gz \
    -o ${subject}_brain_t1w_N4.nii.gz \
    -v 0
  printf " \n ${subject} N4 Correction Complete \n "

  # Brain extraction
  antsBrainExtraction.sh \
    -d 3 \
    -k 1 \
    -z 0 \
    -c 3x1x2x3 \
    -a ${subject}_brain_t1w_N4.nii.gz \
    -e ${oasisPath}/T_template0.nii.gz \
    -m ${oasisPath}/T_template0_BrainCerebellumProbabilityMask.nii.gz \
    -f ${oasisPath}/T_template0_BrainCerebellumRegistrationMask.nii.gz \
    -o ${inputPath}/${subject}/brain/${subject}_brain_t1w_N4

  # Clean up extra output
  rm *Warp.* *Affine* *Tmp.* *0.*

  # # Change the name to something pithier
  # mv ${inputPath}/${subject}/3DT1/${subject}_3DT1_N4BrainExtractionBrain.nii.gz ${inputPath}/${subject}/3DT1/${subject}_3DT1_N4_Brain.nii.gz

  printf " \n Creating ${subject} Brain Extraction Quality Control Images "
  # create mask
  ThresholdImage 3 ${subject}_brain_t1w_N4BrainExtractionSegmentation.nii.gz segmentationMask.nii.gz 0 0 0 1
  # create RGB from segmentation
  # ConvertScalarImageToRGB 3 ${subject}_3DT1_N4BrainExtractionSegmentation.nii.gz segmentationRgb.nii.gz none custom ${qcPath}/snapColormap.txt 0 6
  ConvertScalarImageToRGB 3 ${subject}_brain_t1w_N4BrainExtractionSegmentation.nii.gz segmentationRgb.nii.gz none jet 0 6


  # create tiled mosaic in each orientation
  for dim in 0 1 2
  do
    printf " ${dim} \n "
    printf "${qcPath}/${subject}_${dim}.png "
    CreateTiledMosaic -i ${subject}_brain_t1w_N4.nii.gz -r segmentationRgb.nii.gz -o ${qcPath}/${subject}_brain_t1w_${dim}.png -a 0.3 -t -1x-1 -p mask -s [3,mask,mask] -x segmentationMask.nii.gz -d ${dim}
  done

  # stop the timer
  timer_stop="$(date +"Date : %d/%m/%Y Time : %H.%M.%S")"
  printf " \n \n ${subject} t1w Preparation Complete \n Started: \n ${timer_start} \n Finished: \n ${timer_stop} \n \n "

done
##################################################################################### BRAIN 3DT1 #####################################################################################
# ECOMM


# BCOMM
##################################################################################### BRAIN MWI #####################################################################################
for subject in ${subjects}
do

  # Change into the subject folder
  cd ${inputPath}/${subject}/brain

  # for each mwi dataset
  for mwi in $(ls  *brain_calipr_1.nii* *brain_calipr_2.nii* *brain_grase.nii* *brain_calipr_gematch.* *brain_calipr_ge.* )
  do
    printf " \n Beginning ${subject} ${mwi%.*} Preparation \n \n \n "
    # echo ${mwi%.*}

    # grab E1
    fslroi \
      ${mwi} \
      "${mwi/.nii/_E1.nii}" \
      0 1

    # N4 correction
    N4BiasFieldCorrection \
      -d 3 \
      -i ${mwi%.*}_E1.nii.gz \
      -o ${mwi%.*}_E1_N4.nii.gz

    # Take echo1 to power of 2 to better replicate T1 weighting
    fslmaths \
      ${mwi%.*}_E1_N4.nii.gz \
      -sqr ${mwi%.*}_E1_N4_T1rep.nii.gz

    # Note that ${mwi%.*} gives filename without extension

    # USE -c 3x1x2x3 \ for T1 weighted
    # USE -c 3x3x2x1 \ for T2 weighted
    # Brain Extract MWI 
    antsBrainExtraction.sh \
      -d 3 \
      -k 0 \
      -z 0 \
      -c 3x1x2x3 \
      -a ${mwi%.*}_E1_N4_T1rep.nii.gz \
      -e ${oasisPath}/T_template0.nii.gz \
      -m ${oasisPath}/T_template0_BrainCerebellumProbabilityMask.nii.gz \
      -f ${oasisPath}/T_template0_BrainCerebellumRegistrationMask.nii.gz \
      -o ./${mwi%.*}_E1_N4_T1rep

    # # Clean up extra output
    # rm *Warp.* *Affine* *Tmp.* *0.*

    printf " \n Creating ${mwi%.*} Brain Extraction Quality Control Images "

    # create mask
    ThresholdImage 3 ${mwi%.*}_E1_N4_T1repBrainExtractionMask.nii.gz segmentationMask.nii.gz 0 0 0 1
    # create RGB from segmentation
    ConvertScalarImageToRGB 3 ${mwi%.*}_E1_N4_T1repBrainExtractionMask.nii.gz segmentationRgb.nii.gz none custom ${qcPath}/snapColormap.txt 0 6 

    # create tiled mosaic in each orientation
    for dim in 0 1 2
    do
      CreateTiledMosaic \
        -i ${mwi%.*}_E1_N4_T1rep.nii.gz \
        -r segmentationRgb.nii.gz \
        -o ${qcPath}/${mwi%.*}_${dim}.png \
        -a 0.3 \
        -t -1x-1 \
        -p mask \
        -s [1,mask,mask] \
        -x segmentationMask.nii.gz \
        -d ${dim}
    done
  done
  printf " \n \n ${mwi%.*} Preparation Complete \n \n "

done
##################################################################################### BRAIN MWI #####################################################################################
# ECOMM


# BCOMM
################################################################################ Register MWI -> 3DT1 ################################################################################
for subject in ${subjects}
do
  # Change into the subject folder
  cd ${inputPath}/${subject}/brain

  # for each mwi dataset
  for mwi in $(ls  *brain_calipr_1.nii* *brain_calipr_2.nii* *brain_grase.nii* *brain_calipr_gematch.* *brain_calipr_ge.* )
  do

    # start the timer
    timer_start="$(date +"Date : %d/%m/%Y Time : %H.%M.%S")"
    printf " \n Beginning ${subject} ${mwi%.*} <-> 3DT1 Registration: \n ${timer_start} \n "

    # Register
    antsRegistrationSyN.sh \
      -d 3 \
      -f ${subject}_brain_t1w_N4BrainExtractionBrain.nii.gz \
      -m ${mwi%.*}_E1_N4_T1repBrainExtractionBrain.nii.gz \
      -t a \
      -z 1 \
      -j 0 \
      -p d \
      -n ${cores} \
      -x ${subject}_brain_t1w_N4BrainExtractionMask.nii.gz \
      -o ${mwi%.*}_E1_N4_T1repBrain

    # apply sharper 3DT1 mask to warped GRASE
    fslmaths \
      ${mwi%.*}_E1_N4_T1repBrainWarped.nii.gz \
      -mas ${subject}_brain_t1w_N4BrainExtractionMask.nii.gz \
      ${mwi%.*}_E1_N4_T1repBrainWarpedMasked.nii.gz

    printf " \n Creating ${subject} ${mwi%.*} <-> 3DT1 Registration Quality Control Images "

    # # create mask
    ThresholdImage 3 ${subject}_brain_t1w_N4BrainExtractionMask.nii.gz segmentationMask_wholebrain.nii.gz 0 0 0 1

    for dim in 0 1 2
    do
      # create 3DT1 tiled mosaic
      CreateTiledMosaic -i ${subject}_brain_t1w_N4BrainExtractionBrain.nii.gz -r segmentationRgb.nii.gz -a 0.0 -o ${qcPath}/Reg_${mwi%.*}to3DT1_${dim}_3DT1.png -d ${dim} -p mask -s 5 -x segmentationMask_wholebrain.nii.gz

      # create GRASE tiled mosaic
      CreateTiledMosaic -i ${mwi%.*}_E1_N4_T1repBrainWarpedMasked.nii.gz -r segmentationRgb.nii.gz -a 0.0 -o ${qcPath}/Reg_${mwi%.*}to3DT1_${dim}_MWI.png -d ${dim} -p mask -s 5 -x segmentationMask_wholebrain.nii.gz
    done

    # stop the timer
    timer_stop="$(date +"Date : %d/%m/%Y Time : %H.%M.%S")"
    printf " \n \n ${subject} ${mwi%.*} <-> 3DT1 Registration Complete \n Started: \n ${timer_start} \n Finished: \n ${timer_stop} \n \n "
  done
done
################################################################################ Register MWI -> 3DT1 ################################################################################
# ECOMM


# BCOMM
######################################################################## Warp MWI Images+Maps to 3DT1 Space ########################################################################
for subject in ${subjects}
do
  # Change into the subject folder
  cd ${inputPath}/${subject}/brain

  # for each mwi dataset
  for mwi in $(ls  *brain_calipr_1.nii* *brain_calipr_2.nii* *brain_grase.nii* *brain_calipr_gematch.* *brain_calipr_ge.* )
  do

    printf " \n Beginning ${mwi%.*} Images+Maps warping to 3DT1 \n "
    
    ########## Warp full volume, and create N4 correcting version
    # Redo N4 correction, saving bias field
    N4BiasFieldCorrection \
      -d 3 \
      -i ${mwi%.*}_E1.nii.gz \
      -o [ ${mwi%.*}_E1_N4.nii.gz, ${mwi%.*}_E1_N4biasfield.nii.gz ]

    # make N4 corrected version
    fslmaths \
      ${mwi} \
      -div ${mwi%.*}_E1_N4biasfield.nii.gz \
      ${mwi%.*}_N4.nii.gz
    # warp
    antsApplyTransforms \
      -d 3 \
      -e 3 \
      -i ${mwi%.*}_N4.nii.gz \
      -r ${subject}_brain_t1w_N4BrainExtractionBrain.nii.gz \
      -t ${mwi%.*}_E1_N4_T1repBrain0GenericAffine.mat \
      -o ${mwi%.*}_N4_t1wWarped.nii.gz

    ##### images/masks
    towarp=' E1 E1_N4 E1_N4_T1rep E1_N4_T1repBrainExtractionMask '
    for image in ${towarp}
    do
      antsApplyTransforms \
        -d 3 \
        -i ${mwi%.*}_${image}.nii.gz \
        -r ${subject}_brain_t1w_N4BrainExtractionBrain.nii.gz \
        -t ${mwi%.*}_E1_N4_T1repBrain0GenericAffine.mat \
        -o ${mwi%.*}_${image}_t1wWarped.nii.gz

        # mask with sharper 3DT1 brain mask
        fslmaths \
          ${mwi%.*}_${image}_t1wWarped.nii.gz \
          -mas ${subject}_brain_t1w_N4BrainExtractionMask.nii.gz \
          ${mwi%.*}_${image}_t1wWarped.nii.gz

        # masked version, within WMandGM mask
        fslmaths \
          ${mwi%.*}_${image}_t1wWarped.nii.gz \
          -mas ${subject}_brain_t1w_MaskWMandGM.nii.gz \
          ${mwi%.*}_${image}_t1wWarpedMasked.nii.gz
    done

    ##### quantitative maps
    for iter in {0..2}
    do
      towarp=' MWF IET2 ALPHA '
      for image in ${towarp}
      do
        antsApplyTransforms \
          -d 3 \
          -i ${mwi%.*}_DK${iter}_${image}.nii \
          -r ${subject}_brain_t1w_N4BrainExtractionBrain.nii.gz \
          -t ${mwi%.*}_E1_N4_T1repBrain0GenericAffine.mat \
          -o ${mwi%.*}_DK${iter}_${image}_t1wWarped.nii.gz

        # mask with sharper 3DT1 brain mask
        fslmaths \
          ${mwi%.*}_DK${iter}_${image}_t1wWarped.nii.gz \
          -mas ${subject}_brain_t1w_N4BrainExtractionMask.nii.gz \
          ${mwi%.*}_DK${iter}_${image}_t1wWarped.nii.gz

        # masked version, within WMandGM mask
        fslmaths \
          ${mwi%.*}_DK${iter}_${image}_t1wWarped.nii.gz \
          -mas ${subject}_brain_t1w_MaskWMandGM.nii.gz \
          ${mwi%.*}_DK${iter}_${image}_t1wWarpedMasked.nii.gz
      done
    done

    # Create mask for extracting ROI results from nonzero \
    fslmaths \
      ${mwi%.*}_E1_N4_T1repBrainWarped.nii.gz \
      -mas ${subject}_brain_t1w_N4BrainExtractionMask.nii.gz \
      -bin ${mwi%.*}_E1_N4_T1repBrainWarpedMasked_Bin.nii.gz

    # Create non-zero MWF mask to check ROI size
    fslmaths \
      ${mwi%.*}_DK0_MWF_t1wWarpedMasked.nii.gz \
      -bin ${mwi%.*}_DK0_MWF_t1wWarpedMasked_Bin.nii.gz

    printf " \n Completed ${mwi%.*} Images+Maps warping to 3DT1 \n "
  done
  printf " \n ${subject} qMRI map warping complete \n \n "
done
######################################################################## Warp MWI Images+Maps to 3DT1 Space ########################################################################
# ECOMM


# BCOMM
######################################################################## Create MWI Difference Maps ########################################################################
for subject in ${subjects}
do
  # Change into the subject folder
  cd ${inputPath}/${subject}/brain

  printf " \n ${subject} difference maps \n "

  ##### images (for SNR measurement), loop over E1 and E1_N4
  towarp=' E1 E1_N4 '
  for image in ${towarp}
  do
    # Calculate difference
    fslmaths \
      ${subject}_brain_calipr_1_${image}_t1wWarpedMasked.nii.gz \
      -sub ${subject}_brain_calipr_2_${image}_t1wWarpedMasked.nii.gz \
      -mas ${subject}_brain_t1w_N4BrainExtractionMask.nii.gz \
      ${subject}_brain_calipr_1sub2_${image}_t1wWarpedMasked.nii.gz

    # absolute value version
    fslmaths \
      ${subject}_brain_calipr_1sub2_${image}_t1wWarpedMasked.nii.gz \
      -abs ${subject}_brain_calipr_1sub2abs_${image}_t1wWarpedMasked.nii.gz
  done

  ##### quantitative maps
  for iter in {0..2}
  do
    towarp=' MWF IET2 ALPHA '
    for image in ${towarp}
    do

      # mask with sharper 3DT1 brain mask
      fslmaths \
        ${subject}_brain_calipr_1_DK${iter}_${image}_t1wWarpedMasked.nii.gz \
        -sub ${subject}_brain_calipr_2_DK${iter}_${image}_t1wWarpedMasked.nii.gz \
        -mas ${subject}_brain_t1w_N4BrainExtractionMask.nii.gz \
        ${subject}_brain_calipr_1sub2_DK${iter}_${image}_t1wWarpedMasked.nii.gz
      # absolute value version
      fslmaths \
        ${subject}_brain_calipr_1sub2_DK${iter}_${image}_t1wWarpedMasked.nii.gz \
        -abs ${subject}_brain_calipr_1sub2abs_DK${iter}_${image}_t1wWarpedMasked.nii.gz
    done
  done

done
######################################################################## Create MWI Difference Maps ########################################################################
# ECOMM


# BCOMM
######################################################################## Create C05 GE/GEMATCH MWI Difference Maps ########################################################################
for subject in ${subjects}
do
  # Change into the subject folder
  cd ${inputPath}/${subject}/brain

  printf " \n ${subject} difference maps \n "

  img1='calipr_ge'
  img2='calipr_gematch'


  ##### quantitative maps
  for iter in {0..2}
  do
    towarp=' MWF IET2 ALPHA '
    for image in ${towarp}
    do
      # only masked to brain
      fslmaths \
        ${subject}_brain_${img1}_DK${iter}_${image}_t1wWarped.nii.gz \
        -sub ${subject}_brain_${img2}_DK${iter}_${image}_t1wWarped.nii.gz \
        -mas ${subject}_brain_t1w_N4BrainExtractionMask.nii.gz \
        ${subject}_brain_${img1}sub${img2}_DK${iter}_${image}_t1wWarped.nii.gz
      # absolute value version
      fslmaths \
        ${subject}_brain_${img1}sub${img2}_DK${iter}_${image}_t1wWarped.nii.gz \
        -abs ${subject}_brain_${img1}sub${img2}abs_DK${iter}_${image}_t1wWarped.nii.gz

      # masked to WMandGM
      fslmaths \
        ${subject}_brain_${img1}_DK${iter}_${image}_t1wWarpedMasked.nii.gz \
        -sub ${subject}_brain_${img2}_DK${iter}_${image}_t1wWarpedMasked.nii.gz \
        -mas calipr_repro_c05_brain_t1w_MaskWMandGM.nii.gz \
        ${subject}_brain_${img1}sub${img2}_DK${iter}_${image}_t1wWarpedMasked.nii.gz
      # absolute value version
      fslmaths \
        ${subject}_brain_${img1}sub${img2}_DK${iter}_${image}_t1wWarped.nii.gz \
        -abs ${subject}_brain_${img1}sub${img2}abs_DK${iter}_${image}_t1wWarpedMasked.nii.gz

    done
  done

done
######################################################################## Create C05 GE/GEMATCH MWI Difference Maps ########################################################################
# ECOMM




# BCOMM
############################################################################## 3DT1 <-> Template Registrations ##############################################################################
for subject in ${subjects}
do

  printf " \n Beginning ${subject} ROI registrations \n "

  # Change into the subject folder
  cd ${inputPath}/${subject}/brain

  ################### OLD (but still include, for consistency): Ensure tissue masks lie within brain mask, create a WM+GM mask
  fslmaths \
    ${subject}_brain_t1w_N4BrainExtractionWM.nii.gz \
    -mas ${subject}_brain_t1w_N4BrainExtractionMask.nii.gz -bin \
    ${subject}_brain_t1w_N4BrainExtractionWM.nii.gz

  fslmaths \
    ${subject}_brain_t1w_N4BrainExtractionGM.nii.gz \
    -mas ${subject}_brain_t1w_N4BrainExtractionMask.nii.gz -bin \
    ${subject}_brain_t1w_N4BrainExtractionGM.nii.gz

  # make a WM/GM mask
  fslmaths \
    ${subject}_brain_t1w_N4BrainExtractionWM.nii.gz \
    -add ${subject}_brain_t1w_N4BrainExtractionGM.nii.gz -bin \
    ${subject}_brain_t1w_N4BrainExtractionWMandGM.nii.gz

  
  ################### Register to grase48 template
  antsRegistrationSyN.sh \
    -d 3 \
    -f ${subject}_brain_t1w_N4BrainExtractionBrain.nii.gz \
    -m ${grase48Path}/T_template0.nii.gz \
    -x ${subject}_brain_t1w_N4BrainExtractionMask.nii.gz \
    -t s \
    -z 1 \
    -j 0 \
    -p d \
    -n ${cores} \
    -o ${subject}_T_template

  printf " \n Creating ${subject} Quality Control Images "

  # can skip some redundant steps done earlier

  for dim in 1
  do
    # create 3DT1 tiled mosaic
    CreateTiledMosaic \
      -i ${subject}_brain_t1w_N4BrainExtractionBrain.nii.gz \
      -r segmentationRgb.nii.gz \
      -a 0.0 \
      -o ${qcPath}/Temp_3DT1_Reg_${subject}_${dim}_3DT1.png \
      -t -1x-1 \
      -d ${dim} \
      -p mask \
      -s [4,mask+39,mask] \
      -x segmentationMask_wholebrain.nii.gz

    # create tiled mosaic
    CreateTiledMosaic \
      -i ${subject}_T_templateWarped.nii.gz \
      -r segmentationRgb.nii.gz \
      -a 0.0 \
      -o ${qcPath}/Temp_3DT1_Reg_${subject}_${dim}_MNI.png \
      -t -1x-1 \
      -d ${dim} \
      -p mask \
      -s [4,mask+39,mask] \
      -x segmentationMask_wholebrain.nii.gz
  done

  printf " \n ${subject} ROI registrations complete \n \n "
done
############################################################################## 3DT1 <-> Template Registrations ##############################################################################
# ECOMM


# BCOMM
############################################################################## Mask and ROI Prep ##############################################################################
for subject in ${subjects}
do
  
  # Change into the subject folder
  cd ${inputPath}/${subject}/brain


  ############################################################ Prep Masks

  ################### MWI Overlap Masks
  fsladd \
    ${subject}_brain_t1w_MaskMWIoverlap.nii.gz -m \
    *_E1_N4_T1repBrainExtractionMask_t1wWarped.nii.gz
  
  fslmaths \
    ${subject}_brain_t1w_MaskMWIoverlap.nii.gz -thr 1 \
    ${subject}_brain_t1w_MaskMWIoverlap.nii.gz


  ################### UPDATED Tissue Masks: ensuring mask avoid non-brain signal (eroded) and masked to MWI overlap
  fslmaths \
    ${subject}_brain_t1w_N4BrainExtractionMask.nii.gz \
    -kernel 3D -ero -ero \
    -mas ${subject}_brain_t1w_MaskMWIoverlap.nii.gz -bin \
    ${subject}_brain_t1w_MaskBrain.nii.gz

  fslmaths \
    ${subject}_brain_t1w_MaskBrain.nii.gz \
    -mas ${subject}_brain_t1w_N4BrainExtractionWM.nii.gz -bin \
    ${subject}_brain_t1w_MaskWM.nii.gz

  fslmaths \
    ${subject}_brain_t1w_MaskBrain.nii.gz \
    -mas ${subject}_brain_t1w_N4BrainExtractionGM.nii.gz -bin \
    ${subject}_brain_t1w_MaskGM.nii.gz

  fslmaths \
    ${subject}_brain_t1w_MaskWM.nii.gz \
    -add ${subject}_brain_t1w_MaskGM.nii.gz -bin \
    ${subject}_brain_t1w_MaskWMandGM.nii.gz


  ############################################################ Now for the paper ROIs

  cd ${inputPath}/${subject}/brain

  mkdir ${inputPath}/${subject}/brain/ROIs

  rois=' WMandGM WM GM JHU CC Genu Splenium PIC Caudate Putamen Thalamus Cortical_GM Frontal_WM Occipital_WM Parietal_WM Temporal_WM '

  for roi in ${rois}
  do
    printf "Processing ${subject} ROI ${roi} \n"
    ################### Warp ROI to subject 3DT1
    antsApplyTransforms \
      -d 3 \
      -i ${grase48Path}/Manuscript_ROIs/Brain_ROI_${roi}.nii.gz \
      -r ${subject}_brain_t1w_N4BrainExtractionBrain.nii.gz \
      -t ${subject}_T_template*1Warp.nii.gz \
      -t ${subject}_T_template*0GenericAffine.mat \
      -n MultiLabel \
      -o ROIs/ROI_${roi}_Warped.nii.gz

    # apply sharper masks to labels
    fslmaths \
      ROIs/ROI_${roi}_Warped.nii.gz \
      -mas ${subject}_brain_t1w_MaskWMandGM.nii.gz \
      ROIs/ROI_${roi}_Warped_WMandGMMasked.nii.gz

    fslmaths \
      ROIs/ROI_${roi}_Warped.nii.gz \
      -mas ${subject}_brain_t1w_MaskWM.nii.gz \
      ROIs/ROI_${roi}_Warped_WMMasked.nii.gz

    fslmaths \
      ROIs/ROI_${roi}_Warped.nii.gz \
      -mas ${subject}_brain_t1w_MaskGM.nii.gz \
      ROIs/ROI_${roi}_Warped_GMMasked.nii.gz

  done
  printf " \n ${subject} ROI prep complete \n "

done
############################################################################## Mask and ROI Prep ##############################################################################
# ECOMM


############################################################################## EXTRACT ROI Values ##############################################################################


# BCOMM
##################################################################################### Extract for entire ROI
# rois=' WMandGM WM GM JHU CC Genu Splenium PIC Caudate Putamen Thalamus Cortical_GM Frontal_WM Occipital_WM Parietal_WM Temporal_WM '

rois_WMandGM=' WMandGM GM PIC Caudate Putamen Thalamus '
rois_WM=' WM JHU CC Genu Splenium Frontal_WM Occipital_WM Parietal_WM Temporal_WM '
rois_GM=' Cortical_GM '

# rois_WMandGM=' WMandGM Thalamus '
# rois_WM=' WM JHU '
# rois_GM=' Cortical_GM '


rm ${inputPath}/ROIresults.txt
# Create Header
printf "Subject  Data  Metric  ROI  Mask  ROIMask  Mean  Median  Stddev  5perc  95perc  Volume  Volume_mm3  DKiter  \n" >> ${inputPath}/ROIresults.txt

for subject in ${subjects}
do
  for iter in {0..2}
  do
    # Change into the subject folder
    cd ${inputPath}/${subject}/brain

    # for each MWI dataset
    for mwi in $(ls  *brain_calipr_1.nii* *brain_calipr_2.nii* *brain_grase.nii* *brain_calipr_gematch.* *brain_calipr_ge.* )
    do
      imgs=' MWF IET2 '
      # imgs=' MWF '
      for img in ${imgs}
      do
        
        printf " \n Extracting ${mwi%.*}   DK${iter}   ${img} \n \n "
        ###############################   WMandGM
        for roi in ${rois_WMandGM}
        do
          # loop over roi
          printf "${subject}  ${mwi%.*}  ${img}  ${roi}  WMandGM  MWIoverlap  " >> ${inputPath}/ROIresults.txt
          # mean
          fslstats ${mwi%.*}_DK${iter}_${img}_t1wWarpedMasked.nii.gz -k ROIs/ROI_${roi}_Warped_WMandGMMasked.nii.gz -m | tr '\n' ' ' >> ${inputPath}/ROIresults.txt
          # median
          fslstats ${mwi%.*}_DK${iter}_${img}_t1wWarpedMasked.nii.gz -k ROIs/ROI_${roi}_Warped_WMandGMMasked.nii.gz -p 50 | tr '\n' ' ' >> ${inputPath}/ROIresults.txt
          # std
          fslstats ${mwi%.*}_DK${iter}_${img}_t1wWarpedMasked.nii.gz -k ROIs/ROI_${roi}_Warped_WMandGMMasked.nii.gz -s | tr '\n' ' ' >> ${inputPath}/ROIresults.txt
          # 5th percential
          fslstats ${mwi%.*}_DK${iter}_${img}_t1wWarpedMasked.nii.gz -k ROIs/ROI_${roi}_Warped_WMandGMMasked.nii.gz -p 5 | tr '\n' ' ' >> ${inputPath}/ROIresults.txt
          # 95th percential
          fslstats ${mwi%.*}_DK${iter}_${img}_t1wWarpedMasked.nii.gz -k ROIs/ROI_${roi}_Warped_WMandGMMasked.nii.gz -p 95 | tr '\n' ' ' >> ${inputPath}/ROIresults.txt
          # size of ROI in this image
          fslstats ROIs/ROI_${roi}_Warped_WMandGMMasked.nii.gz -V | tr '\n' ' ' >> ${inputPath}/ROIresults.txt
          # iterations of DK analysis
          printf "${iter}  \n" >> ${inputPath}/ROIresults.txt
        done
        ###############################   WM
        for roi in ${rois_WM}
        do
          # loop over roi
          printf "${subject}  ${mwi%.*}  ${img}  ${roi}  WM  MWIoverlap  " >> ${inputPath}/ROIresults.txt
          # mean
          fslstats ${mwi%.*}_DK${iter}_${img}_t1wWarpedMasked.nii.gz -k ROIs/ROI_${roi}_Warped_WMMasked.nii.gz -m | tr '\n' ' ' >> ${inputPath}/ROIresults.txt
          # median
          fslstats ${mwi%.*}_DK${iter}_${img}_t1wWarpedMasked.nii.gz -k ROIs/ROI_${roi}_Warped_WMMasked.nii.gz -p 50 | tr '\n' ' ' >> ${inputPath}/ROIresults.txt
          # std
          fslstats ${mwi%.*}_DK${iter}_${img}_t1wWarpedMasked.nii.gz -k ROIs/ROI_${roi}_Warped_WMMasked.nii.gz -s | tr '\n' ' ' >> ${inputPath}/ROIresults.txt
          # 5th percential
          fslstats ${mwi%.*}_DK${iter}_${img}_t1wWarpedMasked.nii.gz -k ROIs/ROI_${roi}_Warped_WMMasked.nii.gz -p 5 | tr '\n' ' ' >> ${inputPath}/ROIresults.txt
          # 95th percential
          fslstats ${mwi%.*}_DK${iter}_${img}_t1wWarpedMasked.nii.gz -k ROIs/ROI_${roi}_Warped_WMMasked.nii.gz -p 95 | tr '\n' ' ' >> ${inputPath}/ROIresults.txt
          # size of ROI in this image
          fslstats ROIs/ROI_${roi}_Warped_WMMasked.nii.gz -V | tr '\n' ' ' >> ${inputPath}/ROIresults.txt
          # iterations of DK analysis
          printf "${iter}  \n" >> ${inputPath}/ROIresults.txt
        done
        ###############################   GM
        for roi in ${rois_GM}
        do
          # loop over roi
          printf "${subject}  ${mwi%.*}  ${img}  ${roi}  GM  MWIoverlap  " >> ${inputPath}/ROIresults.txt
          # mean
          fslstats ${mwi%.*}_DK${iter}_${img}_t1wWarpedMasked.nii.gz -k ROIs/ROI_${roi}_Warped_GMMasked.nii.gz -m | tr '\n' ' ' >> ${inputPath}/ROIresults.txt
          # median
          fslstats ${mwi%.*}_DK${iter}_${img}_t1wWarpedMasked.nii.gz -k ROIs/ROI_${roi}_Warped_GMMasked.nii.gz -p 50 | tr '\n' ' ' >> ${inputPath}/ROIresults.txt
          # std
          fslstats ${mwi%.*}_DK${iter}_${img}_t1wWarpedMasked.nii.gz -k ROIs/ROI_${roi}_Warped_GMMasked.nii.gz -s | tr '\n' ' ' >> ${inputPath}/ROIresults.txt
          # 5th percential
          fslstats ${mwi%.*}_DK${iter}_${img}_t1wWarpedMasked.nii.gz -k ROIs/ROI_${roi}_Warped_GMMasked.nii.gz -p 5 | tr '\n' ' ' >> ${inputPath}/ROIresults.txt
          # 95th percential
          fslstats ${mwi%.*}_DK${iter}_${img}_t1wWarpedMasked.nii.gz -k ROIs/ROI_${roi}_Warped_GMMasked.nii.gz -p 95 | tr '\n' ' ' >> ${inputPath}/ROIresults.txt
          # size of ROI in this image
          fslstats ROIs/ROI_${roi}_Warped_GMMasked.nii.gz -V | tr '\n' ' ' >> ${inputPath}/ROIresults.txt
          # iterations of DK analysis
          printf "${iter}  \n" >> ${inputPath}/ROIresults.txt
        done
      done
    done
  done
done

#####################################################################################
# ECOMM


# BCOMM
##################################################################################### B1
rois_WM=' WM CC Genu Splenium '
# rois_WM=' CC '

rm ${inputPath}/ROIresults.txt
# Create Header
printf "Subject  Data  Metric  ROI  Mask  ROIMask  Mean  Median  Stddev  5perc  95perc  Volume  Volume_mm3  DKiter  \n" >> ${inputPath}/ROIresults.txt

for subject in ${subjects}
do

  # Change into the subject folder
  cd ${inputPath}/${subject}/brain

  # DKiter ALPHA
  for iter in {0..2}
  do
    # Change into the subject folder
    cd ${inputPath}/${subject}/brain

    # for each MWI dataset
    for mwi in $(ls  *brain_calipr_1.nii* *brain_calipr_2.nii* *brain_grase.nii* *brain_calipr_gematch.* *brain_calipr_ge.* )
    do
      imgs=' ALPHA '
      for img in ${imgs}
      do
        printf " \n Extracting ${mwi%.*}   DK${iter}   ${img} \n \n "
        ###############################   WM
        for roi in ${rois_WM}
        do
          # loop over roi
          printf "${subject}  ${mwi%.*}  ${img}  ${roi}  WM  MWIoverlap  " >> ${inputPath}/ROIresults.txt
          # mean
          fslstats ${mwi%.*}_DK${iter}_${img}_t1wWarpedMasked.nii.gz -k ROIs/ROI_${roi}_Warped_WMMasked.nii.gz -m | tr '\n' ' ' >> ${inputPath}/ROIresults.txt
          # median
          fslstats ${mwi%.*}_DK${iter}_${img}_t1wWarpedMasked.nii.gz -k ROIs/ROI_${roi}_Warped_WMMasked.nii.gz -p 50 | tr '\n' ' ' >> ${inputPath}/ROIresults.txt
          # std
          fslstats ${mwi%.*}_DK${iter}_${img}_t1wWarpedMasked.nii.gz -k ROIs/ROI_${roi}_Warped_WMMasked.nii.gz -s | tr '\n' ' ' >> ${inputPath}/ROIresults.txt
          # 5th percential
          fslstats ${mwi%.*}_DK${iter}_${img}_t1wWarpedMasked.nii.gz -k ROIs/ROI_${roi}_Warped_WMMasked.nii.gz -p 5 | tr '\n' ' ' >> ${inputPath}/ROIresults.txt
          # 95th percential
          fslstats ${mwi%.*}_DK${iter}_${img}_t1wWarpedMasked.nii.gz -k ROIs/ROI_${roi}_Warped_WMMasked.nii.gz -p 95 | tr '\n' ' ' >> ${inputPath}/ROIresults.txt
          # size of ROI in this image
          fslstats ROIs/ROI_${roi}_Warped_WMMasked.nii.gz -V | tr '\n' ' ' >> ${inputPath}/ROIresults.txt
          # iterations of DK analysis
          printf "${iter}  \n" >> ${inputPath}/ROIresults.txt
        done
      done
    done
  done
done
#####################################################################################
# ECOMM


############################################################################## EXTRACT ROI Values ##############################################################################











