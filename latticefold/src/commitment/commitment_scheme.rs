use cyclotomic_rings::rings::SuitableRing;
#[cfg(feature = "parallel")]
use rayon::prelude::*;
use stark_rings::{
    balanced_decomposition::decompose_balanced_vec,
    cyclotomic_ring::{CRT, ICRT},
    OverField,
};

use super::homomorphic_commitment::Commitment;
use crate::{
    ark_base::*, commitment::CommitmentError, decomposition_parameters::DecompositionParams,
};

/// A concrete instantiation of the Ajtai commitment scheme.
/// Contains a random Ajtai matrix for the corresponding Ajtai parameters
/// `C` is the length of commitment vectors or, equivalently, the number of rows of the Ajtai matrix.
/// `W` is the length of witness vectors or, equivalently, the number of columns of the Ajtai matrix.
/// `NTT` is a suitable cyclotomic ring.
#[derive(Clone, Debug)]
pub struct AjtaiCommitmentScheme<const C: usize, NTT: OverField> {
    matrix: Vec<Vec<NTT>>,
    W: usize,
}

impl<const C: usize, NTT: OverField> TryFrom<(Vec<Vec<NTT>>, usize)>
    for AjtaiCommitmentScheme<C, NTT>
{
    type Error = CommitmentError;

    fn try_from(value: (Vec<Vec<NTT>>, usize)) -> Result<Self, Self::Error> {
        let (matrix, w) = value;

        if matrix.len() != C || matrix[0].len() != w {
            return Err(CommitmentError::WrongAjtaiMatrixDimensions(
                matrix.len(),
                matrix[0].len(),
                C,
                w,
            ));
        }

        let mut ajtai_matrix: Vec<Vec<NTT>> = Vec::with_capacity(C);

        for row in matrix.into_iter() {
            let len = row.len();

            if len != w {
                return Err(CommitmentError::WrongAjtaiMatrixDimensions(C, len, C, w));
            }
            ajtai_matrix.push(row)
        }

        Ok(Self {
            matrix: ajtai_matrix,
            W: w,
        })
    }
}

impl<const C: usize, NTT: OverField> AjtaiCommitmentScheme<C, NTT> {
    /// Returns a random Ajtai commitment matrix
    pub fn rand<Rng: rand::Rng + ?Sized>(rng: &mut Rng, w: usize) -> Self {
        Self {
            matrix: vec![vec![NTT::rand(rng); w]; C],
            W: w
        }
    }
}

impl<const C: usize, NTT: SuitableRing> AjtaiCommitmentScheme<C, NTT> {
    /// Commit to a witness in the NTT form.
    /// The most basic one just multiplies by the matrix.
    pub fn commit_ntt(&self, f: &[NTT], w: usize) -> Result<Commitment<C, NTT>, CommitmentError> {
        if f.len() != w {
            return Err(CommitmentError::WrongWitnessLength(f.len(), w));
        }

        let commitment: Vec<NTT> = cfg_iter!(self.matrix)
            .map(|row| {
                row.iter()
                    .zip(f.iter())
                    .fold(NTT::zero(), |acc, (row_j, f_j)| acc + *row_j * f_j)
            })
            .collect();

        Ok(Commitment::from_vec_raw(commitment))
    }

    /// Commit to a witness in the coefficient form.
    /// Performs NTT on each component of the witness and then does Ajtai commitment.
    pub fn commit_coeff<P: DecompositionParams>(
        &self,
        f: Vec<NTT::CoefficientRepresentation>,
        w: usize
    ) -> Result<Commitment<C, NTT>, CommitmentError> {
        if f.len() != w {
            return Err(CommitmentError::WrongWitnessLength(f.len(), w));
        }

        self.commit_ntt(&CRT::elementwise_crt(f), w)
    }

    /// Takes a coefficient form witness, decomposes it vertically in radix-B,
    /// i.e. computes a preimage G_B^{-1}(w), and Ajtai commits to the result.
    pub fn decompose_and_commit_coeff<P: DecompositionParams>(
        &self,
        f: &[NTT::CoefficientRepresentation],
        w: usize
    ) -> Result<Commitment<C, NTT>, CommitmentError> {
        let f = decompose_balanced_vec(f, P::B, P::L)
            .into_iter()
            .flatten()
            .collect::<Vec<_>>();

        self.commit_coeff::<P>(f, w)
    }

    /// Takes an NTT form witness, transforms it into the coefficient form,
    /// decomposes it vertically in radix-B, i.e.
    /// computes a preimage G_B^{-1}(w), and Ajtai commits to the result.
    pub fn decompose_and_commit_ntt<P: DecompositionParams>(
        &self,
        wit: Vec<NTT>,
        w: usize
    ) -> Result<Commitment<C, NTT>, CommitmentError> {
        let coeff: Vec<NTT::CoefficientRepresentation> = ICRT::elementwise_icrt(wit);

        self.decompose_and_commit_coeff::<P>(&coeff, w)
    }
}

#[cfg(test)]
mod tests {
    use cyclotomic_rings::rings::GoldilocksRingNTT;
    use stark_rings::OverField;

    use super::{AjtaiCommitmentScheme, CommitmentError};
    use crate::ark_base::*;

    pub(crate) fn generate_ajtai<const C: usize, NTT: OverField>(w: usize
    ) -> Result<AjtaiCommitmentScheme<C, NTT>, CommitmentError> {
        let mut matrix = Vec::<Vec<NTT>>::new();

        for i in 0..C {
            let mut row = Vec::<NTT>::new();
            for j in 0..w {
                row.push(NTT::from((i * w + j) as u128));
            }
            matrix.push(row)
        }

        AjtaiCommitmentScheme::try_from((matrix, w))
    }

    #[test]
    fn test_commit_ntt() -> Result<(), CommitmentError> {
        const WITNESS_SIZE: usize = 1 << 15;
        const OUTPUT_SIZE: usize = 9;

        let ajtai_data: AjtaiCommitmentScheme<OUTPUT_SIZE, GoldilocksRingNTT> =
            generate_ajtai(WITNESS_SIZE)?;
        let witness: Vec<_> = (0..(1 << 15)).map(|_| 2_u128.into()).collect();

        let committed = ajtai_data.commit_ntt(&witness, WITNESS_SIZE)?;

        for (i, &x) in committed.as_ref().iter().enumerate() {
            let expected: u128 =
                ((WITNESS_SIZE) * (2 * i * WITNESS_SIZE + (WITNESS_SIZE - 1))) as u128;
            assert_eq!(x, expected.into());
        }

        Ok(())
    }
}
