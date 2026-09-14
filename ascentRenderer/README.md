- Replace case.udf into your udf file.
- The udf file search for `c*/*0.f*`
    - Put all your checkpoint files into folder `c1`
    - For example, the fld files will be `c1/case0.f%5d`
- This do all the rendering in the UDF_Setup and exit right after
- Add rendering pipeline to `ascent.yaml`

You will still need the usual .par and .re2 files as it still needs the mesh and par file for the full setup.

In actual runs, the rendering call can be instead called in UDF_ExecuteStep instead for in-situ rendering.

### On Aurora
See example `s.bin`
```
module use /soft/modulefiles/
module load ascent/release/v0.9.5
export NEKRS_ASCENT_INSTALL_DIR="/soft/visualization/ascent/release/v0.9.5/ascent-checkout/"
```
- Some operations are broken on Aurora, [issue](https://github.com/Nek5000/nekRS_HPCsupport/pull/11)
    - Working: contour, slice
    - Broken: clip, clip_with_field


