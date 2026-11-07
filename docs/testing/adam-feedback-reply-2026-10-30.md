Subject: PyMaCRUISE feedback — geometry review and next diagnostic

Hi Adam,

Thanks for the feedback archive. I confirmed that this is the intended 0.1.6
package and that the inputs match our local test case. The final volumes passed
the original-T1-grid checks, with integer segmentation labels.

For the original reference available here, the affine difference is rounding:
the largest corner displacement is about 0.000015 mm. That cannot explain the
visible superior-brain disagreement. Original raw MaCRUISE VTKs use voxel
indices, while Python's final surfaces use physical RAS millimetres, so those
files need a coordinate conversion for a direct mesh-viewer comparison. The
numeric comparison already performs that conversion and still finds anatomical
differences, including a central mean surface distance of about 0.434 mm.

There is also a convergence issue: your first TOADS pass reached 50 iterations
without the energy threshold, whereas our local run converged at 19. The
archive correctly flags that check as failed. I'm investigating the numerical
divergence; the underlying cause is still unconfirmed.

I've prepared a small review/diagnostic tools ZIP: [INSERT DOWNLOAD LINK]

Please extract it and activate the same 0.1.6 environment. From the extracted
tools folder, run the following with your existing output locations:

python make-review-outputs.py --run "/path/to/outputs" --original-output "/path/to/MaCRUISEv3.2.0/Output" --output-dir "/path/to/new-review-folder"
python collect-feedback.py --run "/path/to/outputs"

This uses the existing results and does not rerun reconstruction. In the new
review folder, compare original_SEG_refine.nii.gz with python_SEG_refine.nii.gz
over the input T1 using nearest-neighbor label display. Compare the paired
original_*_RAS_mm.vtk and python_*_RAS_mm.vtk surfaces. These copies use common
geometry and retain the existing labels and mesh faces.

Please return the refreshed outputs-feedback.zip. It adds aggregate image
fingerprints and geometry details, with no images or mesh arrays. Please also
send stages/n4/attempt-01/t1-n4.nii.gz from this same run so I can compare the
first upstream image directly. Keeping that run in place is sufficient; there
is no need for another full run at this point.

Thanks,
[Your name]
