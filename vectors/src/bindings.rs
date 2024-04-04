use numpy::{PyArray2, PyReadonlyArray1, PyReadonlyArray2, ToPyArray};
use pyo3::prelude::*;

#[pymodule]
#[allow(non_snake_case)]
pub fn vectors<'py>(_py: Python<'py>, m: &'py PyModule) -> PyResult<()> {
    #[pyfn(m)]
    fn createVectors_Matrix<'py>(
        py: Python<'py>,
        PS: PyReadonlyArray2<f64>,
        T: f64,
        l: f64,
        Umax: f64,
        offset: f64,
        inpoFact: u64,
        initialVector: PyReadonlyArray1<f64>,
    ) -> PyResult<&'py PyArray2<f64>> {
        Ok(crate::create_vectors(
            PS.as_array(),
            T,
            l,
            Umax,
            offset,
            inpoFact,
            initialVector.to_owned_array(),
        )
        .to_pyarray(py))
    }

    Ok(())
}
