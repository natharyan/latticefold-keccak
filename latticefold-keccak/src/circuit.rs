use std::{fmt::Debug, time::Instant};

use ark_bls12_381::Fr;
use ark_relations::r1cs::{self, ConstraintSystemRef, Field};
use ark_serialize::{CanonicalSerialize, Compress};
use ark_std::{rand::Rng, vec::Vec, UniformRand};
use arkworks_keccak::{
    constraints::{KeccakCircuit, KeccakMode},
    util::{bytes_to_bitvec, keccak256, sha3_256, shake_128, shake_256},
};
use cyclotomic_rings::{
    challenge_set::LatticefoldChallengeSet,
    rings::{BabyBearChallengeSet, BabyBearRingNTT, SuitableRing},
};
use latticefold::{
    arith::{Arith, Witness, CCCS, CCS, LCCCS},
    commitment::AjtaiCommitmentScheme,
    decomposition_parameters::DecompositionParams,
    nifs::{
        linearization::{LFLinearizationProver, LinearizationProver},
        NIFSProver, NIFSVerifier,
    },
    transcript::{self, poseidon::PoseidonTranscript},
};

use crate::util::{ConstraintSystemExt, field_vec_to_ring_vec};

pub fn z_split<F: Field>(cs: ConstraintSystemRef<F>) -> (usize, usize, Vec<F>, Vec<F>) {
    let public_inputs = cs.ret_instance();
    let witnesses = cs.ret_witness();
    assert_eq!(public_inputs[0], F::one());
    println!("Public inputs: {} + 1", public_inputs.len() - 1);
    println!("Witnesses: {}", witnesses.len());
    (
        public_inputs.len() - 1,
        witnesses.len(),
        public_inputs[1..].to_vec(),
        witnesses,
    )
}

fn ret_ccs<
    const C: usize, // rows
    const W: usize, // columns
    P: DecompositionParams,
    R: Clone + UniformRand + Debug + SuitableRing,
    F: Field,
>(
    x_r1cs: Vec<F>,
    w_r1cs: Vec<F>,
    wit_len: usize,
    r1cs_rows: usize,
) -> (
    CCCS<C, R>,
    Witness<R>,
    CCS<R>,
    AjtaiCommitmentScheme<C, W, R>,
) {
    let mut rng = ark_std::test_rng();
    let new_r1cs_rows = if P::L == 1 && (wit_len > 0 && (wit_len & (wit_len - 1)) == 0) {
        r1cs_rows - 2
    } else {
        r1cs_rows
    };
    let x_ccs: Vec<R> = field_vec_to_ring_vec(&x_r1cs); // TODO: field to ring basefield
    let w_ccs: Vec<R> = field_vec_to_ring_vec(&w_r1cs); // TODO: field to ring basefield
    let one = R::one();

    let mut z = vec![one];
    z.extend(&x_ccs);
    z.extend(&w_ccs);
    let ccs: CCS<R> = r1cs_to_ccs::<R, x_r1cs.len(), wit_len, W>(&z, P::L, new_r1cs_rows); // TODO: r1cs_to_ccs
    ccs.check_relation(&z).expect("R1CS invalid!");

    let scheme: AjtaiCommitmentScheme<C, W, R> = AjtaiCommitmentScheme::rand(&mut rng);
    let wit: Witness<R> = Witness::from_w_ccs::<P>(w_ccs);

    let cm_i: CCCS<C, R> = CCCS {
        cm: wit.commit::<C, W, P>(&scheme).unwrap(),
        x_ccs,
    };

    (cm_i, wit, ccs, scheme)
}

fn setup_environment<
    const C: usize,
    RqNTT: SuitableRing,
    DP: DecompositionParams,
    const W: usize,
    CS: LatticefoldChallengeSet<RqNTT>,
    F: Field,
>(
    cs: ConstraintSystemRef<F>,
) -> (
    LCCCS<C, RqNTT>,
    Witness<RqNTT>,
    CCCS<C, RqNTT>,
    Witness<RqNTT>,
    CCS<RqNTT>,
    AjtaiCommitmentScheme<C, W, RqNTT>,
) {
    let (x_len, wit_len, x_r1cs, w_r1cs) = z_split(cs);
    let r1cs_rows = x_len + wit_len + 1;
    let (cm_i, wit, ccs, scheme) = ret_ccs::<C, W, DP, RqNTT, F>(x_r1cs, w_r1cs.clone(), wit_len, r1cs_rows);
    let ring_w_css: Vec<RqNTT> = w_r1cs.clone()
        .into_iter()
        .map(|bit| {
            if bit == F::one() {
                RqNTT::from(1u64)
            } else {
                RqNTT::from(0u64)
            }
        })
        .collect();
    let wit_acc = Witness::from_w_ccs::<DP>(ring_w_css);

    let mut transcript = PoseidonTranscript::<RqNTT, CS>::default();

    let (acc, _) = LFLinearizationProver::<_, PoseidonTranscript<RqNTT, CS>>::prove(
        &cm_i,
        &wit_acc,
        &mut transcript,
        &ccs,
    )
    .expect("Failed to generate linearization proof");

    (acc, wit_acc, cm_i, wit, ccs, scheme)
}
