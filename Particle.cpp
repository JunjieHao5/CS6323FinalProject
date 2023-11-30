#include "Particle.h"

Particle::Particle(){}

Particle::Particle(const Vector2f& pos, const Vector2f& vel, float mass, float lame_lambda, float lame_mu){
	position.setData(pos);
	velocity.setData(vel);
	this->mass = mass;
	lambda = lame_lambda;
	mu = lame_mu;
	//To start out with, we assume the deformation gradient is zero
	//Or in other words, all particle velocities are the same
	def_elastic.loadIdentity();
	def_plastic.loadIdentity();
	svd_e.setData(1, 1);
	svd_w.loadIdentity();
	svd_v.loadIdentity();
	polar_r.loadIdentity();
	polar_s.loadIdentity();
}

Particle::~Particle(){}

void Particle::updatePos(){
	//Simple euler integration
	position += TIMESTEP*velocity;
}


void Particle::updateGradient(){
	//So, initially we make all updates elastic
	velocity_gradient *= TIMESTEP;
	velocity_gradient.diag_sum(1);
	def_elastic.setData(velocity_gradient * def_elastic);
}


void Particle::applyPlasticity(){
	Matrix2f f_all = def_elastic * def_plastic;

	//We compute the SVD decomposition
	//The singular values (basically a scale transform) tell us if 
	//the particle has exceeded critical stretch/compression
	def_elastic.svd(&svd_w, &svd_e, &svd_v);

	Matrix2f svd_v_trans = svd_v.transpose();

	//Clamp singular values to within elastic region
	for (int i=0; i<2; i++){
		if (svd_e[i] < CRIT_COMPRESS)
			svd_e[i] = CRIT_COMPRESS;
		else if (svd_e[i] > CRIT_STRETCH)
			svd_e[i] = CRIT_STRETCH;
	}

#if ENABLE_IMPLICIT
	//Compute polar decomposition, from clamped SVD
	polar_r.setData(svd_w*svd_v_trans);
	polar_s.setData(svd_v);
	polar_s.diag_product(svd_e);
	polar_s.setData(polar_s*svd_v_trans);
#endif
	
	//Recompute elastic and plastic gradient
	//We're basically just putting the SVD back together again
	Matrix2f v_cpy(svd_v), w_cpy(svd_w);
	v_cpy.diag_product_inv(svd_e);
	w_cpy.diag_product(svd_e);
	def_plastic = v_cpy*svd_w.transpose()*f_all;
	def_elastic = w_cpy*svd_v.transpose();
}


const Matrix2f Particle::energyDerivative(){
	//Adjust lame parameters to account for hardening
	float harden = exp(HARDENING*(1-def_plastic.determinant())),
		Je = svd_e.product();
	//This is the co-rotational term
	Matrix2f temp = 2*mu*(def_elastic - svd_w*svd_v.transpose())*def_elastic.transpose();
	//Add in the primary contour term
	temp.diag_sum(lambda*Je*(Je-1));
	//Add hardening and volume
	return volume * harden * temp;
}


#if ENABLE_IMPLICIT
const Vector2f Particle::deltaForce(const Vector2f& u, const Vector2f& weight_grad) {
    // Calculate delta(Fe) for the elastic deformation gradient
    Matrix2f del_elastic = TIMESTEP * u.outer_product(weight_grad) * def_elastic;
    
    // Skip calculations if delta(Fe) is negligible
    if (del_elastic.isNegligible(MATRIX_EPSILON))
        return Vector2f(0);

    // Compute skew symmetric part of R^T*dF - dF^TR
    float y = (polar_r[0][0] * del_elastic[1][0] + polar_r[1][0] * del_elastic[1][1]) -
              (polar_r[0][1] * del_elastic[0][0] + polar_r[1][1] * del_elastic[0][1]);

    // Solve for x in the linear system derived from MS + SM
    float x = y / (polar_s[0][0] + polar_s[1][1]);

    // Compute deltaR = R*(R^T*dR)
    Matrix2f del_rotate(
        -polar_r[1][0] * x, polar_r[0][0] * x,
        -polar_r[1][1] * x, polar_r[0][1] * x
    );
    
    // Compute cofactor matrix of F, JF^-T
    Matrix2f cofactor = def_elastic.cofactor();
        
    // Compute delta(JF^-T) as the cofactor of delta(F)
    Matrix2f del_cofactor = del_elastic.cofactor();

    // Calculate "A" for the force computation
    Matrix2f Ap = del_elastic - del_rotate;
    Ap *= 2 * mu;
    cofactor *= cofactor.frobeniusInnerProduct(del_elastic);
    del_cofactor *= (def_elastic.determinant() - 1);
    cofactor += del_cofactor;
    cofactor *= lambda;
    Ap += cofactor;
    
    // Combine all components for the final force calculation
    return volume * (Ap * (def_elastic.transpose() * weight_grad));
}
#endif

