use numpy::{PyArray2, PyReadonlyArray1, PyReadonlyArray2, ToPyArray};
use pyo3::prelude::*;

#[pymodule]
#[allow(non_snake_case)]
pub fn totalRotMatrix<'py>(_py: Python<'py>, Rtot: &'py PyModule) -> PyResult<()> {
    #[pyfn(Rtot)]
    fn create_Rtot<'py>(
        py: Python<'py>,
        PS: PyReadonlyArray2<f64>,
        T: f64,
        Umax: f64,
        offset: f64,
        inpoFact: u64,
    ) -> PyResult<&'py PyArray2<f64>> {
        Ok(crate::create_Rtot(
            PS.as_array(),
            T,
            Umax,
            offset,
            inpoFact,
        )
        .to_pyarray(py))
    }

    Ok(())
}
