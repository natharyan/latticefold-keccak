use ark_ff::PrimeField;
use ark_relations::r1cs::{ConstraintSystemRef, Field};
use cyclotomic_rings::rings::SuitableRing;
use num_bigint::BigInt;

pub trait ConstraintSystemExt<F: Field> {
    fn ret_instance(&self) -> Vec<F>;
    fn ret_witness(&self) -> Vec<F>;
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
}


// TODO: field to ring basefield
pub fn field_vec_to_ring_vec<R,F>(elems: &[F]) -> Vec<R>
where 
R: SuitableRing, 
F: Field,
R::BaseRing: From<F>,
{
    elems
        .iter()
        .map(|f| {
            let mut coeffs = vec![R::BaseRing::ZERO; R::dimension()];
            coeffs[0] = R::BaseRing::from(*f);
            R::from(coeffs)
        })
        .collect()
}
