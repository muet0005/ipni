# AUTO-GENERATED by tools/checkspecs.py - DO NOT EDIT
from nipype.testing import assert_equal
from nipype.interfaces.slicer.segmentation.specialized import BRAINSROIAuto

def test_BRAINSROIAuto_inputs():
    input_map = dict(ROIAutoDilateSize=dict(argstr='--ROIAutoDilateSize %f',
    ),
    args=dict(argstr='%s',
    ),
    closingSize=dict(argstr='--closingSize %f',
    ),
    environ=dict(nohash=True,
    usedefault=True,
    ),
    ignore_exception=dict(nohash=True,
    usedefault=True,
    ),
    inputVolume=dict(argstr='--inputVolume %s',
    ),
    numberOfThreads=dict(argstr='--numberOfThreads %d',
    ),
    otsuPercentileThreshold=dict(argstr='--otsuPercentileThreshold %f',
    ),
    outputClippedVolumeROI=dict(argstr='--outputClippedVolumeROI %s',
    hash_files=False,
    ),
    outputROIMaskVolume=dict(argstr='--outputROIMaskVolume %s',
    hash_files=False,
    ),
    outputVolumePixelType=dict(argstr='--outputVolumePixelType %s',
    ),
    terminal_output=dict(mandatory=True,
    nohash=True,
    ),
    thresholdCorrectionFactor=dict(argstr='--thresholdCorrectionFactor %f',
    ),
    )
    inputs = BRAINSROIAuto.input_spec()

    for key, metadata in input_map.items():
        for metakey, value in metadata.items():
            yield assert_equal, getattr(inputs.traits()[key], metakey), value

def test_BRAINSROIAuto_outputs():
    output_map = dict(outputClippedVolumeROI=dict(),
    outputROIMaskVolume=dict(),
    )
    outputs = BRAINSROIAuto.output_spec()

    for key, metadata in output_map.items():
        for metakey, value in metadata.items():
            yield assert_equal, getattr(outputs.traits()[key], metakey), value
