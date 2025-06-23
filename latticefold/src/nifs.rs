//! The NIFS module defines the behaviour of the [LatticeFold](https://eprint.iacr.org/2024/257.pdf) protocol
//!
//! NIFS = Non Interactive Folding Scheme

use ark_ff::{Field, PrimeField};
use ark_serialize::{CanonicalDeserialize, CanonicalSerialize};
use ark_std::{marker::PhantomData, vec::Vec};
use cyclotomic_rings::rings::SuitableRing;
use stark_rings::OverField;

use self::{decomposition::*, error::LatticefoldError, folding::*, linearization::*};
use crate::{
    arith::{error::CSError, Witness, CCCS, CCS, LCCCS},
    commitment::AjtaiCommitmentScheme,
    decomposition_parameters::DecompositionParams,
    transcript::{Transcript, TranscriptWithShortChallenges},
};

pub mod decomposition;
pub mod error;
pub mod folding;
pub mod linearization;

#[cfg(test)]
mod tests;

/// `C` is the length of Ajtai commitment vectors.
/// `NTT` is a cyclotomic ring in the NTT form.
#[derive(Clone, CanonicalSerialize, CanonicalDeserialize)]
pub struct LFProof<const C: usize, NTT: OverField> {
    pub linearization_proof: LinearizationProof<NTT>,
    pub decomposition_proof_l: DecompositionProof<C, NTT>,
    pub decomposition_proof_r: DecompositionProof<C, NTT>,
    pub folding_proof: FoldingProof<NTT>,
}

/// `C` is the length of commitment vectors or, equivalently, the number of rows of the Ajtai matrix.
/// `W` is the length of witness vectors or, equivalently, the number of columns of the Ajtai matrix.
/// `NTT` is a suitable cyclotomic ring.
/// `P` is the decomposition parameters.
/// `T` is the FS-transform transcript.
pub struct NIFSProver<const C: usize, NTT, P, T> {
    _r: PhantomData<NTT>,
    _p: PhantomData<P>,
    _t: PhantomData<T>,
}

impl<
        const C: usize,
        NTT: SuitableRing,
        P: DecompositionParams,
        T: TranscriptWithShortChallenges<NTT>,
    > NIFSProver<C, NTT, P, T>
{
    pub fn prove(
        acc: &LCCCS<C, NTT>,
        w_acc: &Witness<NTT>,
        cm_i: &CCCS<C, NTT>,
        w_i: &Witness<NTT>,
        transcript: &mut impl TranscriptWithShortChallenges<NTT>,
        ccs: &CCS<NTT>,
        scheme: &AjtaiCommitmentScheme<C, NTT>,
        w: usize
    ) -> Result<(LCCCS<C, NTT>, Witness<NTT>, LFProof<C, NTT>), LatticefoldError<NTT>> {
        sanity_check::<NTT, P>(ccs)?;

        absorb_public_input::<NTT, C>(acc, cm_i, transcript);

        let (linearized_cm_i, linearization_proof) =
            LFLinearizationProver::<_, T>::prove(cm_i, w_i, transcript, ccs)?;
        let (mz_mles_l, decomposed_lcccs_l, decomposed_wit_l, decomposition_proof_l) =
            LFDecompositionProver::<_, T>::prove::<C, P>(acc, w_acc, transcript, ccs, scheme, w)?;
        let (mz_mles_r, decomposed_lcccs_r, decomposed_wit_r, decomposition_proof_r) =
            LFDecompositionProver::<_, T>::prove::<C, P>(
                &linearized_cm_i,
                w_i,
                transcript,
                ccs,
                scheme,
                w
            )?;

        let (mz_mles, lcccs, wit_s) = {
            let mut lcccs = decomposed_lcccs_l;
            let mut lcccs_r = decomposed_lcccs_r;
            lcccs.append(&mut lcccs_r);

            let mut wit_s = decomposed_wit_l;
            let mut wit_s_r = decomposed_wit_r;
            wit_s.append(&mut wit_s_r);

            let mut mz_mles = mz_mles_l;
            let mut mz_mles_r = mz_mles_r;
            mz_mles.append(&mut mz_mles_r);
            (mz_mles, lcccs, wit_s)
        };

        let (folded_lcccs, wit, folding_proof) =
            LFFoldingProver::<_, T>::prove::<C, P>(&lcccs, wit_s, transcript, ccs, &mz_mles)?;

        Ok((
            folded_lcccs,
            wit,
            LFProof {
                linearization_proof,
                decomposition_proof_l,
                decomposition_proof_r,
                folding_proof,
            },
        ))
    }
}

/// `C` is the length of commitment vectors or, equivalently, the number of rows of the Ajtai matrix.
/// `W` is the length of witness vectors or, equivalently, the number of columns of the Ajtai matrix.
/// `NTT` is a suitable cyclotomic ring.
/// `P` is the decomposition parameters.
/// `T` is the FS-transform transcript.
pub struct NIFSVerifier<const C: usize, NTT, P, T> {
    _r: PhantomData<NTT>,
    _p: PhantomData<P>,
    _t: PhantomData<T>,
}

impl<
        const C: usize,
        NTT: SuitableRing,
        P: DecompositionParams,
        T: TranscriptWithShortChallenges<NTT>,
    > NIFSVerifier<C, NTT, P, T>
{
    pub fn verify(
        acc: &LCCCS<C, NTT>,
        cm_i: &CCCS<C, NTT>,
        proof: &LFProof<C, NTT>,
        transcript: &mut impl TranscriptWithShortChallenges<NTT>,
        ccs: &CCS<NTT>,
    ) -> Result<LCCCS<C, NTT>, LatticefoldError<NTT>> {
        sanity_check::<NTT, P>(ccs)?;

        absorb_public_input::<NTT, C>(acc, cm_i, transcript);

        let linearized_cm_i = LFLinearizationVerifier::<_, T>::verify(
            cm_i,
            &proof.linearization_proof,
            transcript,
            ccs,
        )?;
        let decomposed_acc = LFDecompositionVerifier::<_, T>::verify::<C, P>(
            acc,
            &proof.decomposition_proof_l,
            transcript,
            ccs,
        )?;
        let decomposed_cm_i = LFDecompositionVerifier::<_, T>::verify::<C, P>(
            &linearized_cm_i,
            &proof.decomposition_proof_r,
            transcript,
            ccs,
        )?;

        let lcccs_s = {
            let mut decomposed_acc = decomposed_acc;
            let mut decomposed_cm_i = decomposed_cm_i;

            decomposed_acc.append(&mut decomposed_cm_i);

            decomposed_acc
        };

        Ok(LFFoldingVerifier::<NTT, T>::verify::<C, P>(
            &lcccs_s,
            &proof.folding_proof,
            transcript,
            ccs,
        )?)
    }
}

fn sanity_check<NTT: SuitableRing, DP: DecompositionParams>(
    ccs: &CCS<NTT>,
) -> Result<(), LatticefoldError<NTT>> {
    if ccs.m != usize::max((ccs.n - ccs.l - 1) * DP::L, ccs.m).next_power_of_two() {
        return Err(CSError::InvalidSizeBounds(ccs.m, ccs.n, DP::L).into());
    }

    Ok(())
}

fn absorb_public_input<NTT: SuitableRing, const C: usize>(
    acc: &LCCCS<C, NTT>,
    cm_i: &CCCS<C, NTT>,
    transcript: &mut impl Transcript<NTT>,
) {
    transcript.absorb_field_element(&<NTT::BaseRing as Field>::from_base_prime_field(
        <NTT::BaseRing as Field>::BasePrimeField::from_be_bytes_mod_order(b"acc"),
    ));

    transcript.absorb_slice(&acc.r);
    transcript.absorb_slice(&acc.v);
    transcript.absorb_slice(acc.cm.as_ref());
    transcript.absorb_slice(&acc.u);
    transcript.absorb_slice(&acc.x_w);
    transcript.absorb(&acc.h);

    transcript.absorb_field_element(&<NTT::BaseRing as Field>::from_base_prime_field(
        <NTT::BaseRing as Field>::BasePrimeField::from_be_bytes_mod_order(b"cm_i"),
    ));

    transcript.absorb_slice(cm_i.cm.as_ref());
    transcript.absorb_slice(&cm_i.x_ccs);
}
