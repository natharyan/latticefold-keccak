use ark_std::test_rng;
use cyclotomic_rings::{
    challenge_set::LatticefoldChallengeSet,
    rings::{BabyBearChallengeSet, FrogChallengeSet, GoldilocksChallengeSet, StarkChallengeSet},
};
use ark_bls12_381::Fr;
use num_traits::One;
use rand::Rng;
use stark_rings::cyclotomic_ring::models::{
    babybear::RqNTT as BabyBearRqNTT, frog_ring::RqNTT as FrogRqNTT,
    goldilocks::RqNTT as GoldilocksRqNTT, stark_prime::RqNTT as StarkRqNTT,
};

use super::*;
use crate::{
    arith::{r1cs::get_test_z_split, tests::get_test_ccs, utils::mat_vec_mul},
    commitment::AjtaiCommitmentScheme,
    decomposition_parameters::{
        test_params::{BabyBearDP, FrogDP, GoldilocksDP, StarkDP},
        DecompositionParams,
    },
    nifs::linearization::utils::{sumcheck_polynomial_comb_fn, SqueezeBeta},
    transcript::poseidon::PoseidonTranscript,
};

const C: usize = 4;
const WIT_LEN: usize = 4;
fn setup_test_environment<
    RqNTT: SuitableRing,
    DP: DecompositionParams,
    const C: usize,
    const W: usize,
>(
    input: Option<usize>,
) -> (
    Witness<RqNTT>,
    CCCS<C, RqNTT>,
    CCS<RqNTT>,
    AjtaiCommitmentScheme<C, RqNTT>,
) {
    let ccs = get_test_ccs::<RqNTT>(W, DP::L);
    let mut rng = test_rng();
    let (_, x_ccs, w_ccs) = get_test_z_split::<RqNTT>(input.unwrap_or(rng.gen_range(0..64)));
    let scheme = AjtaiCommitmentScheme::rand(&mut rng, W);

    let wit = Witness::from_w_ccs::<DP>(w_ccs);
    let cm_i = CCCS {
        cm: wit.commit::<C, DP>(&scheme, WIT_LEN).unwrap(),
        x_ccs,
    };

    (wit, cm_i, ccs, scheme)
}

#[test]
fn test_compute_z_ccs() {
    type RqNTT = StarkRqNTT;
    type DP = StarkDP;
    const W: usize = WIT_LEN * DP::L;
    let (wit, cm_i, _, scheme) = setup_test_environment::<RqNTT, DP, C, W>(None);

    let z_ccs = <CCCS<C, RqNTT> as Instance<RqNTT, Fr>>::get_z_vector(&cm_i,&wit.w_ccs);

    // Check z_ccs structure
    assert_eq!(z_ccs.len(), cm_i.x_ccs.len() + 1 + wit.w_ccs.len());
    assert_eq!(z_ccs[cm_i.x_ccs.len()], RqNTT::one());

    // Check commitment
    assert_eq!(cm_i.cm, wit.commit::<C, StarkDP>(&scheme, W).unwrap());
}

#[test]
fn test_construct_polynomial() {
    type RqNTT = GoldilocksRqNTT;
    type CS = GoldilocksChallengeSet;
    type DP = GoldilocksDP;
    const W: usize = WIT_LEN * DP::L;
    let (wit, cm_i, ccs, _) = setup_test_environment::<RqNTT, DP, C, W>(None);

    let z_ccs = <CCCS<C, RqNTT> as Instance<RqNTT, Fr>>::get_z_vector(&cm_i,&wit.w_ccs);

    let mut transcript = PoseidonTranscript::<RqNTT, CS>::default();
    let (_, g_degree, mz_mles) =
        LFLinearizationProver::<RqNTT, PoseidonTranscript<RqNTT, CS>>::construct_polynomial_g(
            &z_ccs,
            &mut transcript,
            &ccs,
        )
        .unwrap();

    // Check dimensions
    assert_eq!(mz_mles.len(), ccs.t);

    // Check degree of g
    assert!(g_degree <= ccs.q + 1)
}

#[test]
fn test_generate_sumcheck() {
    type RqNTT = FrogRqNTT;
    type CS = FrogChallengeSet;
    type DP = FrogDP;
    const W: usize = WIT_LEN * DP::L;
    let (wit, cm_i, ccs, _) = setup_test_environment::<RqNTT, DP, C, W>(None);

    let z_ccs = <CCCS<C, RqNTT> as Instance<RqNTT, Fr>>::get_z_vector(&cm_i,&wit.w_ccs);

    let mut transcript = PoseidonTranscript::<RqNTT, CS>::default();
    let (g_mles, g_degree, _) =
        LFLinearizationProver::<RqNTT, PoseidonTranscript<RqNTT, CS>>::construct_polynomial_g(
            &z_ccs,
            &mut transcript,
            &ccs,
        )
        .unwrap();

    let comb_fn = |vals: &[RqNTT]| -> RqNTT { sumcheck_polynomial_comb_fn::<RqNTT>(vals, &ccs) };

    let (_, point_r) =
        LFLinearizationProver::<RqNTT, PoseidonTranscript<RqNTT, CS>>::generate_sumcheck_proof(
            &mut transcript,
            g_mles,
            ccs.s,
            g_degree,
            comb_fn,
        )
        .unwrap();

    // Check dimensions
    assert_eq!(point_r.len(), ccs.s);
}

fn prepare_test_vectors<RqNTT: SuitableRing, CS: LatticefoldChallengeSet<RqNTT>, F: PrimeField>(
    wit: &Witness<RqNTT>,
    cm_i: &CCCS<C, RqNTT>,
    ccs: &CCS<RqNTT>,
) -> (Vec<RqNTT>, Vec<DenseMultilinearExtension<RqNTT>>) {
    let z_ccs = <CCCS<C, RqNTT> as Instance<RqNTT, F>>::get_z_vector(cm_i,&wit.w_ccs);

    let mut transcript = PoseidonTranscript::<RqNTT, CS>::default();
    let (g_mles, g_degree, Mz_mles) =
        LFLinearizationProver::<RqNTT, PoseidonTranscript<RqNTT, CS>>::construct_polynomial_g(
            &z_ccs,
            &mut transcript,
            ccs,
        )
        .unwrap();

    let comb_fn = |vals: &[RqNTT]| -> RqNTT { sumcheck_polynomial_comb_fn::<RqNTT>(vals, ccs) };

    let (_, point_r) =
        LFLinearizationProver::<RqNTT, PoseidonTranscript<RqNTT, CS>>::generate_sumcheck_proof(
            &mut transcript,
            g_mles,
            ccs.s,
            g_degree,
            comb_fn,
        )
        .unwrap();

    (point_r, Mz_mles)
}

#[test]
fn test_compute_v() {
    type RqNTT = BabyBearRqNTT;
    type CS = BabyBearChallengeSet;
    type DP = BabyBearDP;
    const W: usize = WIT_LEN * DP::L;


    // Setup shared test state
    let (wit, cm_i, ccs, _) = setup_test_environment::<RqNTT, DP, C, W>(None);
    let (point_r, Mz_mles) = prepare_test_vectors::<RqNTT, CS, Fr>(&wit, &cm_i, &ccs);

    // Compute actual v vector
    let (point_r, v, _) =
        LFLinearizationProver::<RqNTT, PoseidonTranscript<RqNTT, CS>>::compute_evaluation_vectors(
            &wit, &point_r, &Mz_mles,
        )
        .unwrap();

    // Compute expected v vector (witness evaluations)
    let expected_v =
        evaluate_mles::<RqNTT, _, _, LinearizationError<RqNTT>>(&wit.f_hat, &point_r).unwrap();

    // Validate
    assert_eq!(point_r.len(), ccs.s, "point_r length mismatch");
    assert!(!v.is_empty(), "v vector should not be empty");
    assert_eq!(
        v, expected_v,
        "v vector doesn't match expected witness evaluations"
    );
}

#[test]
fn test_compute_u() {
    type RqNTT = FrogRqNTT;
    type CS = FrogChallengeSet;
    type DP = FrogDP;
    const W: usize = WIT_LEN * DP::L;
    // Setup shared test state
    let (wit, cm_i, ccs, _) = setup_test_environment::<RqNTT, DP, C, W>(None);
    let (point_r, Mz_mles) = prepare_test_vectors::<RqNTT, CS, Fr>(&wit, &cm_i, &ccs);

    // Compute actual u vector
    let (point_r, _, u) =
        LFLinearizationProver::<RqNTT, PoseidonTranscript<RqNTT, CS>>::compute_evaluation_vectors(
            &wit, &point_r, &Mz_mles,
        )
        .unwrap();

    // Compute expected u vector
    let z_ccs = <CCCS<C, RqNTT> as Instance<RqNTT, Fr>>::get_z_vector(&cm_i,&wit.w_ccs);
    let expected_Mz_mles: Vec<DenseMultilinearExtension<RqNTT>> = ccs
        .M
        .iter()
        .map(|M| {
            DenseMultilinearExtension::from_evaluations_vec(ccs.s, mat_vec_mul(M, &z_ccs).unwrap())
        })
        .collect();
    let expected_u = compute_u(&expected_Mz_mles, &point_r).unwrap();

    // Validate
    assert_eq!(point_r.len(), ccs.s, "point_r length mismatch");
    assert!(!u.is_empty(), "u vector should not be empty");
    assert_eq!(u, expected_u, "u vector doesn't match expected evaluations");
}

#[test]
fn test_full_prove() {
    type RqNTT = GoldilocksRqNTT;
    type CS = GoldilocksChallengeSet;
    type DP = GoldilocksDP;
    const W: usize = WIT_LEN * DP::L;
    let (wit, cm_i, ccs, _) = setup_test_environment::<RqNTT, DP, C, W>(None);
    let mut transcript = PoseidonTranscript::<RqNTT, CS>::default();

    let (lcccs, proof) = LFLinearizationProver::<RqNTT, PoseidonTranscript<RqNTT, CS>>::prove::<C,Fr>(
        &cm_i,
        &wit,
        &mut transcript,
        &ccs,
    )
    .unwrap();

    assert_eq!(lcccs.r.len(), ccs.s);
    assert_eq!(lcccs.v.len(), proof.v.len());
    assert_eq!(lcccs.u.len(), proof.u.len());
}

#[test]
fn test_verify_sumcheck_proof() {
    type RqNTT = StarkRqNTT;
    type CS = StarkChallengeSet;
    type DP = StarkDP;
    const W: usize = WIT_LEN * DP::L;
    let (wit, cm_i, ccs, _) = setup_test_environment::<RqNTT, DP, C, W>(None);
    let mut prove_transcript = PoseidonTranscript::<RqNTT, CS>::default();

    // Generate proof
    let (lcccs, proof) = LFLinearizationProver::<RqNTT, PoseidonTranscript<RqNTT, CS>>::prove::<C, Fr>(
        &cm_i,
        &wit,
        &mut prove_transcript,
        &ccs,
    )
    .unwrap();

    // We need to recreate the exact same transcript state
    let mut verify_transcript = PoseidonTranscript::<RqNTT, CS>::default();

    // Generate beta challenges to match prover's transcript state
    let _ = verify_transcript.squeeze_beta_challenges(ccs.s);

    let result =
        LFLinearizationVerifier::<RqNTT, PoseidonTranscript<RqNTT, CS>>::verify_sumcheck_proof(
            &proof,
            &mut verify_transcript,
            &ccs,
        );

    // Instead of unwrapping, handle the result
    match result {
        Ok((point_r, _)) => {
            assert_eq!(point_r.len(), ccs.s);
            // We know that point_r from lcccs is valid
            assert_eq!(point_r, lcccs.r);
        }
        Err(e) => panic!("Sumcheck verification failed: {:?}", e),
    }
}

#[test]
fn test_verify_evaluation_claim() {
    type RqNTT = BabyBearRqNTT;
    type CS = BabyBearChallengeSet;
    type DP = BabyBearDP;
    const W: usize = WIT_LEN * DP::L;
    let (wit, cm_i, ccs, _) = setup_test_environment::<RqNTT, DP, C, W>(None);
    let mut transcript = PoseidonTranscript::<RqNTT, CS>::default();

    // Generate proof
    let (_, proof) = LFLinearizationProver::<RqNTT, PoseidonTranscript<RqNTT, CS>>::prove::<C, Fr>(
        &cm_i,
        &wit,
        &mut transcript,
        &ccs,
    )
    .unwrap();

    // Reset transcript and generate verification data
    let mut transcript = PoseidonTranscript::<RqNTT, CS>::default();
    let beta_s = transcript.squeeze_beta_challenges(ccs.s);

    let (point_r, s) =
        LFLinearizationVerifier::<RqNTT, PoseidonTranscript<RqNTT, CS>>::verify_sumcheck_proof(
            &proof,
            &mut transcript,
            &ccs,
        )
        .unwrap();

    // Test the evaluation claim verification
    let result =
        LFLinearizationVerifier::<RqNTT, PoseidonTranscript<RqNTT, CS>>::verify_evaluation_claim(
            &beta_s, &point_r, s, &proof, &ccs,
        );

    assert!(result.is_ok());
}

#[test]
fn test_prepare_verifier_output() {
    type RqNTT = FrogRqNTT;
    type CS = FrogChallengeSet;
    type DP = FrogDP;
    const W: usize = WIT_LEN * DP::L;
    let (wit, cm_i, ccs, _) = setup_test_environment::<RqNTT, DP, C, W>(None);
    let mut transcript = PoseidonTranscript::<RqNTT, CS>::default();

    let (_, proof) = LFLinearizationProver::<RqNTT, PoseidonTranscript<RqNTT, CS>>::prove::<C, Fr>(
        &cm_i,
        &wit,
        &mut transcript,
        &ccs,
    )
    .unwrap();

    let point_r = vec![RqNTT::one(); ccs.s]; // Example point_r

    let lcccs =
        LFLinearizationVerifier::<RqNTT, PoseidonTranscript<RqNTT, CS>>::prepare_verifier_output::<C, Fr>(
            &cm_i,
            point_r.clone(),
            &proof,
        );

    // Verify final state structure
    assert_eq!(lcccs.r, point_r);
    assert_eq!(lcccs.v, proof.v);
    assert_eq!(lcccs.u, proof.u);
    assert_eq!(lcccs.cm, cm_i.cm);
    assert_eq!(lcccs.x_w, cm_i.x_ccs);
    assert_eq!(lcccs.h, RqNTT::one());
}

#[test]
fn test_verify_invalid_proof() {
    type RqNTT = GoldilocksRqNTT;
    type CS = GoldilocksChallengeSet;
    type DP = GoldilocksDP;
    const W: usize = WIT_LEN * DP::L;
    let (wit, cm_i, ccs, _) = setup_test_environment::<RqNTT, DP, C, W>(None);
    let mut transcript = PoseidonTranscript::<RqNTT, CS>::default();

    let (_, mut proof) = LFLinearizationProver::<RqNTT, PoseidonTranscript<RqNTT, CS>>::prove::<C, Fr>(
        &cm_i,
        &wit,
        &mut transcript,
        &ccs,
    )
    .unwrap();

    // Corrupt the proof
    if !proof.u.is_empty() {
        proof.u[0] += RqNTT::one();
    }

    // Reset transcript for verification
    let mut transcript = PoseidonTranscript::<RqNTT, CS>::default();

    let result = LFLinearizationVerifier::<RqNTT, PoseidonTranscript<RqNTT, CS>>::verify::<C, Fr>(
        &cm_i,
        &proof,
        &mut transcript,
        &ccs,
    );

    assert!(result.is_err());
}
