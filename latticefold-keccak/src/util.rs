use ark_relations::r1cs::{ConstraintSystemRef, Field, ConstraintMatrices};
use cyclotomic_rings::rings::SuitableRing;
use stark_rings::Ring;
use stark_rings_linalg::SparseMatrix;

pub trait ConstraintSystemExt<F: Field> {
    fn ret_instance(&self) -> Vec<F>;
    fn ret_witness(&self) -> Vec<F>;
    fn get_r1cs_matrices<R: Ring + std::convert::From<F>>(&self) -> (SparseMatrix<R>, SparseMatrix<R>, SparseMatrix<R>);
}

impl<F: Field> ConstraintSystemExt<F> for ConstraintSystemRef<F> {
    fn ret_instance(&self) -> Vec<F> {
        let cs = self.borrow().unwrap();
        cs.instance_assignment.clone()
    }

    fn ret_witness(&self) -> Vec<F> {
        let cs = self.borrow().unwrap();
        cs.witness_assignment.clone()
    }

    // TODO: finish this
    fn get_r1cs_matrices<R: Ring + std::convert::From<F>>(&self) -> (SparseMatrix<R>, SparseMatrix<R>, SparseMatrix<R>) {
        // get A, B, C;
        // get matrices for A, B, C from cs; use these to construct sparse matrices.
        let matrices: ConstraintMatrices<F> = self.borrow().unwrap().to_matrices().unwrap();
        // (elem, col_index)
        let a: Vec<Vec<(F, usize)>> = matrices.a;
        let b: Vec<Vec<(F, usize)>> = matrices.b;
        let c: Vec<Vec<(F, usize)>> = matrices.c;
        let a = a.r1csmat_to_sparsematrix();
        let b = b.r1csmat_to_sparsematrix();
        let c = c.r1csmat_to_sparsematrix();
        (a, b, c)
    }
}


pub trait MatrixExt<F: Field, R: Ring> {
    fn r1csmat_to_sparsematrix(&self) -> SparseMatrix<R>;
}

impl<F: Field, R: Ring + std::convert::From<F>> MatrixExt<F, R> for Vec<Vec<(F, usize)>> {
    fn r1csmat_to_sparsematrix(&self) -> SparseMatrix<R> {
        let n_rows = self.len();
        let n_cols = self.iter().flat_map(|row| row.iter().map(|&(_, col)| col)).max().map_or(0, |max_col| max_col + 1);
        let mut matrix = SparseMatrix {
            n_rows: n_rows,
            n_cols: n_cols,
            coeffs: vec![vec![]; n_rows],
        };

        for (i, row) in matrix.coeffs.iter_mut().enumerate() {
            for (elem, j) in &self[i] {
                row.push((field_to_ring(*elem), *j));
            }
        }
        matrix
    }
}

pub fn hadamard<R: Ring>(a: SparseMatrix<R>, b: SparseMatrix<R>, c: SparseMatrix<R>, z: &[R]) -> SparseMatrix<R> {
    unimplemented!()
}

pub fn field_to_ring<R, F>(elem: F) -> R
where
    R: Ring + std::convert::From<F>,
    F: Field,
{
    R::from(elem)
}

pub fn fieldvec_to_ringvec<R, F>(elems: &[F]) -> Vec<R>
where
    R: SuitableRing,
    F: Field,
{
    elems
        .iter()
        .map(|f| {
            let mut coeffs = vec![R::BaseRing::from(0 as u128); R::dimension()];
            // expect boolean witnesses from r1cs-std
            coeffs[0] = if *f == F::zero() {
                R::BaseRing::from(0 as u128)
            } else if *f == F::one() {
                R::BaseRing::from(1 as u128)
            } else {
                panic!("Unexpected field element: {:?}", f);
            };
            R::from(coeffs)
        })
        .collect()
}


