//! Provides operations used for working with constraint systems

#[cfg(feature = "parallel")]
use rayon::iter::{IntoParallelRefIterator, ParallelIterator};
use stark_rings::Ring;
use cyclotomic_rings::rings::SuitableRing;
use stark_rings_linalg::SparseMatrix;

use super::error::CSError as Error;
use crate::ark_base::*;
use ark_ff::PrimeField;

//  Computes the hadamard product of two ring
#[allow(dead_code)]
pub(crate) fn hadamard_vec<R: Ring>(lhs: &[R], rhs: &[R]) -> Vec<R> {
    lhs.iter().zip(rhs).map(|(lhs, rhs)| *lhs * rhs).collect()
}

// Multiplies Vector of rings by another ring
#[allow(dead_code)]
pub(crate) fn vec_value_mul<R: Ring>(lhs: &[R], rhs: &R) -> Vec<R> {
    lhs.iter().map(|lhs_i| *lhs_i * rhs).collect()
}

// Adds two ring vectors
pub(crate) fn vec_add<R: Ring>(a: &[R], b: &[R]) -> Result<Vec<R>, Error> {
    if a.len() != b.len() {
        return Err(Error::LengthsNotEqual(
            "a".to_string(),
            "b".to_string(),
            a.len(),
            b.len(),
        ));
    }
    Ok(a.iter().zip(b.iter()).map(|(x, y)| *x + y).collect())
}

pub(crate) fn vec_scalar_mul<R: Ring>(vec: &[R], c: &R) -> Vec<R> {
    vec.iter().map(|a| *a * c).collect()
}

pub(crate) fn hadamard<R: Ring>(a: &[R], b: &[R]) -> Result<Vec<R>, Error> {
    if a.len() != b.len() {
        return Err(Error::LengthsNotEqual(
            "a".to_string(),
            "b".to_string(),
            a.len(),
            b.len(),
        ));
    }
    Ok(a.iter().zip(b).map(|(a, b)| *a * b).collect())
}

pub(crate) fn mat_vec_mul<R: Ring>(M: &SparseMatrix<R>, z: &[R]) -> Result<Vec<R>, Error> {
    if M.n_cols != z.len() {
        return Err(Error::LengthsNotEqual(
            "M".to_string(),
            "z".to_string(),
            M.n_cols,
            z.len(),
        ));
    }

    Ok(cfg_iter!(M.coeffs)
        .map(|row| row.iter().map(|(value, col_i)| *value * z[*col_i]).sum())
        .collect())
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

pub(crate) fn fieldvec_to_ringvec<R, F>(elems: &[F]) -> Vec<R>
where
    R: SuitableRing,
    F: PrimeField,
{
    elems
        .iter()
        .map(|f| {
            let mut coeffs = vec![R::BaseRing::from(0 as u128); R::dimension()];
            coeffs[0] = field_to_ring(*f);
            R::from(coeffs)
        })
        .collect()
}

#[cfg(test)]
mod tests {
    use ark_ff::Zero;
    use stark_rings::cyclotomic_ring::models::goldilocks::Fq;
    use stark_rings_linalg::sparse_matrix::dense_matrix_to_sparse;

    use super::*;

    #[test]
    fn test_hadamard_vec() {
        let a = [Fq::from(2u64), Fq::from(3u64), Fq::from(4u64)];
        let b = [Fq::from(5u64), Fq::from(6u64), Fq::from(7u64)];
        let result = hadamard_vec(&a, &b);
        let expected = vec![Fq::from(10u64), Fq::from(18u64), Fq::from(28u64)];
        assert_eq!(result, expected);
    }

    // Add similar tests for other functions here...

    #[test]
    fn test_vec_value_mul() {
        let a = [Fq::from(2u64), Fq::from(3u64), Fq::from(4u64)];
        let scalar = Fq::from(2u64);
        let result = vec_value_mul(&a, &scalar);
        let expected = vec![Fq::from(4u64), Fq::from(6u64), Fq::from(8u64)];
        assert_eq!(result, expected);
    }

    #[test]
    fn test_vec_add() {
        let a = [Fq::from(1u64), Fq::from(2u64), Fq::from(3u64)];
        let b = [Fq::from(4u64), Fq::from(5u64), Fq::from(6u64)];
        let result = vec_add(&a, &b);
        let expected = vec![Fq::from(5u64), Fq::from(7u64), Fq::from(9u64)];
        assert_eq!(result.unwrap(), expected);

        // Test error case
        let a = [Fq::from(1u64), Fq::from(2u64)];
        let b = [Fq::from(3u64), Fq::from(4u64), Fq::from(5u64)];
        let result = vec_add(&a, &b);
        assert!(result.is_err());
    }

    #[test]
    fn test_vec_scalar_mul() {
        let vec = [Fq::from(1u64), Fq::from(2u64), Fq::from(3u64)];
        let c = Fq::from(3u64);
        let result = vec_scalar_mul(&vec, &c);
        let expected = vec![Fq::from(3u64), Fq::from(6u64), Fq::from(9u64)];
        assert_eq!(result, expected);
    }

    #[test]
    fn test_hadamard() {
        let a = [Fq::from(2u64), Fq::from(3u64), Fq::from(4u64)];
        let b = [Fq::from(5u64), Fq::from(6u64), Fq::from(7u64)];
        let result = hadamard(&a, &b);
        let expected = vec![Fq::from(10u64), Fq::from(18u64), Fq::from(28u64)];
        assert_eq!(result.unwrap(), expected);

        // Test error case
        let a = [Fq::from(2u64), Fq::from(3u64)];
        let b = [Fq::from(5u64), Fq::from(6u64), Fq::from(7u64)];
        let result = hadamard(&a, &b);
        assert!(result.is_err());
    }

    #[test]
    fn test_mat_vec_mul() {
        // Construct a sparse matrix M
        let dense_matrix = vec![
            vec![Fq::from(1u64), Fq::zero(), Fq::zero()], // Row 0
            vec![Fq::zero(), Fq::from(2u64), Fq::from(1u64)], // Row 1
            vec![Fq::zero(), Fq::zero(), Fq::from(3u64)], // Row 2
        ];

        let M = dense_matrix_to_sparse(dense_matrix);

        let z = [Fq::from(1u64), Fq::from(1u64), Fq::from(1u64)];
        let result = mat_vec_mul(&M, &z);
        let expected = vec![Fq::from(1u64), Fq::from(3u64), Fq::from(3u64)];
        assert_eq!(result.unwrap(), expected);

        // Test error case
        let z = [Fq::from(1u64), Fq::from(1u64)]; // Wrong size vector
        let result = mat_vec_mul(&M, &z);
        assert!(result.is_err());
    }
}
