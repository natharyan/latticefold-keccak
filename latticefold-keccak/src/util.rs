use std::{usize, vec};

use ark_ff::PrimeField;
use ark_relations::r1cs::{ConstraintMatrices, ConstraintSystemRef, Field};
use cyclotomic_rings::rings::SuitableRing;
use stark_rings::Ring;
use stark_rings_linalg::SparseMatrix;

pub trait ConstraintSystemExt<F: PrimeField> {
    fn ret_instance(&self) -> Vec<F>;
    fn ret_witness(&self) -> Vec<F>;
    fn get_r1cs_matrices<R: Ring>(&self) -> (SparseMatrix<R>, SparseMatrix<R>, SparseMatrix<R>);
}

impl<F: PrimeField> ConstraintSystemExt<F> for ConstraintSystemRef<F> {
    fn ret_instance(&self) -> Vec<F> {
        let cs = self.borrow().unwrap();
        cs.instance_assignment.clone()
    }

    fn ret_witness(&self) -> Vec<F> {
        let cs = self.borrow().unwrap();
        cs.witness_assignment.clone()
    }

    fn get_r1cs_matrices<R: Ring>(&self) -> (SparseMatrix<R>, SparseMatrix<R>, SparseMatrix<R>) {
        // get A, B, C;
        // get matrices for A, B, C from cs; use these to construct sparse matrices.
        let matrices: ConstraintMatrices<F> = self.borrow().unwrap().to_matrices().unwrap();
        // (elem, col_index)
        let a: Vec<Vec<(F, usize)>> = matrices.a;
        let b: Vec<Vec<(F, usize)>> = matrices.b;
        let c: Vec<Vec<(F, usize)>> = matrices.c;

        // let a_n_rows = a.len();
        // let a_n_cols = a
        //     .iter()
        //     .flat_map(|row| row.iter().map(|&(_, col)| col))
        //     .max()
        //     .map_or(0, |max_col| max_col + 1);
        // println!("a: ({}, {})", a_n_rows, a_n_cols);

        // let b_n_rows = b.len();
        // let b_n_cols = b
        //     .iter()
        //     .flat_map(|row| row.iter().map(|&(_, col)| col))
        //     .max()
        //     .map_or(0, |max_col| max_col + 1);
        // println!("b: ({}, {})", b_n_rows, b_n_cols);

        // let c_n_rows = c.len();
        let c_n_cols = c
            .iter()
            .flat_map(|row| row.iter().map(|&(_, col)| col))
            .max()
            .map_or(0, |max_col| max_col + 1);
        // println!("c: ({}, {})\n", c_n_rows, c_n_cols);

        let a = a.r1csmat_to_sparsematrix(c_n_cols);
        let b = b.r1csmat_to_sparsematrix(c_n_cols);
        let c = c.r1csmat_to_sparsematrix(c_n_cols);
        (a, b, c)
    }
}

pub trait MatrixExt<F: PrimeField, R: Ring> {
    fn r1csmat_to_sparsematrix(&self, n_cols: usize) -> SparseMatrix<R>;
}

impl<F: PrimeField, R: Ring> MatrixExt<F, R> for Vec<Vec<(F, usize)>> {
    fn r1csmat_to_sparsematrix(&self, n_cols: usize) -> SparseMatrix<R> {
        let n_rows = self.len();
        // assert_eq!(n_rows, n_cols, "non square r1cs matrix");
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

// find one such solution for d such that (d . z) = (a . z)o(b . z)o(c . z) = ((diag(a . z) . diag(b . z)) . c) . z
pub fn hadamard_ret_d<R: Ring>(
    a: SparseMatrix<R>,
    b: SparseMatrix<R>,
    c: SparseMatrix<R>,
    z: &[R],
    rows: usize,
) -> SparseMatrix<R> {
    assert_eq!(
        c.n_cols,
        z.len(),
        "Length of witness vector must be equal to ccs width"
    );
    // matrix width need not match - depends on the witness used in a constraint

    let mut matrix = SparseMatrix {
        n_rows: rows,
        n_cols: z.len(),
        coeffs: vec![vec![]; rows],
    };

    let matvec_mult_a_z = sparse_matrix_mult_vec(&a.coeffs, z);
    let matvec_mult_b_z = sparse_matrix_mult_vec(&b.coeffs, z);
    let mut diag_hadamard_a_b_z = Vec::new();
    for i in 0..rows {
        diag_hadamard_a_b_z.push(matvec_mult_a_z[i] * matvec_mult_b_z[i]);
    }
    for (i, row) in c.coeffs.iter().enumerate() {
        for (elem, col) in row {
            matrix.coeffs[i].push((diag_hadamard_a_b_z[i] * (*elem), *col));
        }
    }

    matrix
}

pub fn field_to_ring<R, F>(elem: F) -> R
where
    R: Ring,
    F: PrimeField,
{
    let bigint: <F as PrimeField>::BigInt = elem.into_bigint();
    let val_bigint: u64 = bigint.as_ref()[0];
    R::from(val_bigint as u128)
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

fn sparse_matrix_mult_vec<R: Ring>(rows: &Vec<Vec<(R, usize)>>, wit_vec: &[R]) -> Vec<R> {
    let mut result = Vec::new();
    for row in rows {
        let mut dotproduct = R::from(0u128);
        for (elem, index) in row {
            dotproduct += *elem * wit_vec[*index];
        }
        result.push(dotproduct);
    }
    result
}

pub fn pad_matrtixrows_to<R: Ring>(
    mut matrix: SparseMatrix<R>,
    num_cols: usize,
) -> SparseMatrix<R> {
    for row in matrix.coeffs.iter_mut() {
        let mut ctr = 0;
        let mut new_elements = Vec::new();
        for (_, col) in row.iter() {
            if *col > ctr {
                new_elements.push((R::from(0 as u128), ctr));
                ctr += 1;
            } else {
                ctr += 1;
            }
        }
        while ctr < num_cols {
            new_elements.push((R::from(0 as u128), ctr));
            ctr += 1;
        }
        row.extend(new_elements);
        row.sort_by(|(_, col1), (_, col2)| col1.cmp(col2));
    }
    matrix
}
