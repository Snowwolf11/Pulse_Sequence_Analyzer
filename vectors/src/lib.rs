use ndarray::{s, Array1, Array2, Array3, ArrayView2, ArrayViewMut2, Axis};

#[cfg(feature = "python")]
mod bindings;

#[allow(non_snake_case)]
pub fn create_vectors(
    PS: ArrayView2<f64>,
    T: f64,
    l: f64,
    Umax: f64,
    offset: f64,
    inpoFact: u64,
    initialVector: Array1<f64>,
) -> Array2<f64> {
    let mut m = Array2::<f64>::zeros((PS.shape()[0] + 1, 3));

    m.slice_mut(s![0, ..]).assign(&(initialVector * l));

    if PS.shape()[0] == 0 {
        return m;
    }

    let angle_factor = -2.0 * std::f64::consts::PI * T / inpoFact as f64;
    let offset_factor = Umax * inpoFact as f64 / 100.0;

    let mut rotation_matrices = Array3::<f64>::zeros((PS.shape()[0], 3, 3));

    let iter = ndarray::Zip::from(rotation_matrices.outer_iter_mut())
        .and(PS.slice(s![.., 0]))
        .and(PS.slice(s![.., 1]));

    #[cfg(feature = "par")]
    {
        iter.par_for_each(|rotation_matrix, &PS0, &PS1| {
            calc_rotation_matix(
                rotation_matrix,
                PS0,
                PS1,
                offset,
                angle_factor,
                offset_factor,
            )
        });
    }
    #[cfg(not(feature = "par"))]
    {
        iter.for_each(|rotation_matrix, &PS0, &PS1| {
            calc_rotation_matix(
                rotation_matrix,
                PS0,
                PS1,
                offset,
                angle_factor,
                offset_factor,
            )
        });
    }

    let vn = m.slice(s![0, ..]).to_owned();
    let mut accum = Array2::<f64>::eye(3);
    let iter = rotation_matrices
        .outer_iter()
        .zip(m.axis_iter_mut(Axis(0)).skip(1));
    for (rotation_matrix, mut m) in iter {
        accum = accum.dot(&rotation_matrix);
        let vn = accum.dot(&vn);
        m.assign(&vn);
    }

    m
}

#[allow(non_snake_case)]
pub fn calc_rotation_matix(
    mut rotation_matrix: ArrayViewMut2<f64>,
    PS0: f64,
    PS1: f64,
    offset: f64,
    angle_factor: f64,
    offset_factor: f64,
) {
    let Ux = PS0 * offset_factor * PS1.to_radians().cos();
    let Uy = PS0 * offset_factor * PS1.to_radians().sin();
    let Uz = offset;

    let norm_n = (Ux.powi(2) + Uy.powi(2) + Uz.powi(2)).sqrt();
    if norm_n != 0.0 {
        let n0 = Ux / norm_n;
        let n1 = Uy / norm_n;
        let n2 = Uz / norm_n;

        let cosa = (angle_factor * norm_n).cos();
        let sina = (angle_factor * norm_n).sin();
        let mcosa = 1.0 - cosa;
        rotation_matrix.assign(&ndarray::arr2(&[
            [
                n0.powi(2) * mcosa + cosa,
                n0 * n1 * mcosa - n2 * sina,
                n0 * n2 * mcosa + n1 * sina,
            ],
            [
                n0 * n1 * mcosa + n2 * sina,
                n1.powi(2) * mcosa + cosa,
                n1 * n2 * mcosa - n0 * sina,
            ],
            [
                n2 * n0 * mcosa - n1 * sina,
                n2 * n1 * mcosa + n0 * sina,
                n2.powi(2) * mcosa + cosa,
            ],
        ]));
    } else {
        rotation_matrix.assign(&Array2::eye(3));
    }
}
