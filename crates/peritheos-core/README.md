# peritheos-core

`peritheos-core` is the dependency-light native numerical core of Peritheos.
It is under active migration and is not yet published as a compatibility-ready
crate.

The crate owns built-in equations, validation, thermodynamic properties, and
inversion. It does not depend on Python, NumPy, PyO3, or a fitting framework.
See the repository's Rust migration contract for the compatibility and
validation gates that apply before release.

Scalar APIs are the reference behavior. Ordered, fail-fast batch evaluation is
available through the dependency-free extension traits:

```rust
use peritheos_core::batch::IsothermalEosBatch;
use peritheos_core::isothermal::BM3;

let eos = BM3::new(10.0, 160.0, 4.0)?;
let pressures = eos.pressures(&[8.0, 9.0, 10.0])?;
# Ok::<(), peritheos_core::EosError>(())
```

The batch API deliberately does not choose a threading or array framework.
Embedders can parallelize at their own boundary; the Python bindings use
thresholded Rayon evaluation for large NumPy arrays.
