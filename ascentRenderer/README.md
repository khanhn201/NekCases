- Replace case.udf with your udf file.
- Add rendering pipeline to `ascent.yaml`

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


