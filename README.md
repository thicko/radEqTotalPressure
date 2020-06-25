# rotEqTotalPressure
Radially varying total pressure boundary condition for OpenFOAM. 
Adapted from `Foam::rotatingTotalPressureFvPatchScalarField`. 
(https://github.com/OpenFOAM/OpenFOAM-dev/tree/master/src/finiteVolume/fields/fvPatchFields/derived/rotatingTotalPressure)

## Installation

1. Clone the repository into `$WM_PROJECT_USER_DIR/src/finiteVolume/fields/fvPatchFields/derived`.
2. Copy the make directory to `$WM_PROJECT_USER_DIR/src/finiteVolumne`.
3. Run `wmake libso` in `$WM_PROJECT_USER_DIR/src/finiteVolume`.

## Use

Tell OpenFOAM to dynamically link the library by including

```
libs
(
    "librotEqTotalPressureFvPatchScalarField.so"
);
```

in `<MY_CASE>/system/controlDict`.

In the file `<MY_CASE>/0/p` or `<MY_CASE>/0/p_rgh`, the boundary condition can be set using e.g.:

```
<PATCH_NAME>
{
    type            rotEqTotalPressure;
    p0              uniform 100000;
    rRef            0.1;
    dp0dr           1000;
    rho             thermo:rho;
    psi             thermo:psi;
    gamma           1.4;
    omega           (0 0 100);
}
```

Property     | Description             | Required    | Default value
:---         |:---                    |:---       |:---
`U`            | velocity field name     | no          | U
`phi`          | flux field name         | no          | phi
`rho`          | density field name      | no          | none
`psi`          | compressibility field name | no       | none
`gamma`        | ratio of specific heats (Cp/Cv) | yes | none
`p0`           | static pressure reference | yes       | none
`omega`        | angular velocity of the frame [rad/s] | yes    | none
`dp0dr`        | radial total pressure gradient | yes | none
`rRef`         | reference radial location for p0 |  yes| none

The boundary condition is currently restricted to rotation around the Z-axis, so `omega` should be specified in the form

```
omega (0 0 w);
```

where `w` is the frame rotational speed in rad/s.


