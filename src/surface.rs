
// ==============================================================================
// PART 1: Traits & Interfaces
// ==============================================================================

pub trait GemSurface {
    fn radius_a(&self) -> f64;
    fn volume(&self) -> f64;
    fn surface_area(&self) -> f64;
    fn parametric_surface(&self, u: f64, v: f64) -> [f64; 3];
    fn implicit_equation(&self, x: f64, y: f64, z: f64) -> f64;
    fn metric_tensor(&self, v: f64) -> (f64, f64) {
        let a = self.radius_a();
        let cos_half_v = (v / 2.0).cos();
        let g_uu = 4.0 * a.powi(2) * cos_half_v.powi(4);
        let g_vv = a.powi(2);
        (g_uu, g_vv)
    }
    fn gaussian_curvature(&self, v: f64) -> f64 {
        let a = self.radius_a();
        let denom = a.powi(2) * (1.0 + v.cos());
        if denom.abs() < 1e-9 { return 0.0; }
        v.cos() / denom
    }
}