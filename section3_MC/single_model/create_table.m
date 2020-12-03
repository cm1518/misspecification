med=median(draws,2);
fifth=prctile(draws,5,2);
ninetyfifth=prctile(draws,95,2);
 
rowLabels = {'$\alpha$', '$\gamma_b$','$\iota_w$','$100\gamma$','$h$','$\lambda_p^{ss}$','$\lambda_w^{ss}$','$L^{ss}$','$\pi^{ss}$','$100 (\beta^{-1}-1)$','$\nu$','$\kappa_p$','$\xi_w$','$\chi$', '$S$','$\phi_{\pi}$','$\phi_{y}$','$\phi_{dy}$','$\rho_R$','$\rho_z$','$\rho_g$','$\rho_{\mu}$','$\rho_{p}$','$\rho_{w}$','$\rho_b$','$\rho_{mp}$','$\theta_p$','$\theta_w$','$\sigma_{mp}$','$\sigma_{z}$','$\sigma_{g}$','$\sigma_{\mu}$','$\sigma_{p}$','$\sigma_{w}$','$\sigma_{b}$'};
   columnLabels = {'5th percentile', 'median','95th percentile'};
  matrix2latex([fifth med ninetyfifth], 'out.tex', 'rowLabels', rowLabels, 'columnLabels', columnLabels, 'alignment', 'c', 'format', '%-6.2f', 'size', 'tiny');