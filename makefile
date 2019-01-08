SHELL=/bin/bash
filename=jpaper
LATEX_ROOT=latex
FIG_ROOT=${LATEX_ROOT}/figures

MATLAB=matlab -nodesktop -nosplash -r
INIT_CMD = addpath functions;
ML_CMD         = "try;                                   \
					${INIT_CMD}                          \
					run $^                               \
						catch E;                         \
						fprintf('%s', E.message);        \
						exit(1);                         \
				  end;                                   \
				  exit(0);"


SCRIPT_ROOT=publish_plotting

fig_FRF=${FIG_ROOT}/G_uz2stage_Eres.svg ${FIG_ROOT}/G_pow_and_current.svg
m_FRF=${SCRIPT_ROOT}/plot_FRF_xdir.m

fig_drift=${FIG_ROOT}/drift_fit.svg
m_drift=${SCRIPT_ROOT}/plot_fit_drift.m

fig_hyst=${FIG_ROOT}/hyst_response_ol.svg
m_hyst=${SCRIPT_ROOT}/plot_hyst_sat_fit.m

fig_gvib=${FIG_ROOT}/pzplot.svg
tbl_gvib=${LATEX_ROOT}/Gvib_data.tex
m_gvib=${SCRIPT_ROOT}/table_pzdata.m

fig_locusCZ=${FIG_ROOT}/lqr_locus_choozezet.svg
m_locusCZ=${SCRIPT_ROOT}/plot_lqr_locus_choosezeta.m

fig_locusCS=${FIG_ROOT}/lqr_locus_constsig_0p9.svg
m_locusCS=${SCRIPT_ROOT}/plot_lqr_locus_constsigma.m


tbl_tsmeans=${LATEX_ROOT}/ts_means_CZ.tex ${LATEX_ROOT}/ts_means_CS.tex
m_tsmeans=${SCRIPT_ROOT}/table_ts_means.m

fig_SBWgmpm=${FIG_ROOT}/GainS_TS_vs_gamma_both.svg \
            ${FIG_ROOT}/BW_TS_vs_gamma_both.svg    \
	        ${FIG_ROOT}/PMGM_vs_gamma_both.svg
m_SBWgmpm=${SCRIPT_ROOT}/plot_gmpm_gainS_both_schemes_newts.m

fig_exp=${FIG_ROOT}/step_exps_min_gam.svg \
	    ${FIG_ROOT}/step_exps_rob_gam.svg \
		${FIG_ROOT}/steps.svg
m_exp=${SCRIPT_ROOT}/plot_experiment_vs_sim.m

assets=${fig_drift}   \
	   ${fig_hyst}    \
	   ${fig_FRF}     \
       ${fig_gvib}    \
       ${tbl_gvib}    \
       ${fig_locusCZ} \
       ${fig_locusCS} \
	   ${fig_SBWgmpm} \
       ${tbl_tsmeans} \
	   ${fig_exp}

# assets=${fig_exp}


pdf:${assets}
	echo "in PDF"
	#echo ${assets}
	cd ${LATEX_ROOT} && make


${fig_drift}: ${m_drift}
	#echo "IN DRIFT"
	${MATLAB} ${ML_CMD}

${fig_hyst}: ${m_hyst}
	${MATLAB} ${ML_CMD}

${fig_FRF}: ${m_FRF}
	${MATLAB} ${ML_CMD}

${fig_gvib} ${tbl_gvib}: ${m_gvib}
	${MATLAB} ${ML_CMD}

${fig_locusCS}: ${m_locusCS}
	${MATLAB} ${ML_CMD}

${fig_locusCZ}: ${m_locusCZ}
	${MATLAB} ${ML_CMD}

${tbl_tsmeans}: ${m_tsmeans}
	${MATLAB} ${ML_CMD}

${fig_SBWgmpm}: ${m_SBWgmpm}
	${MATLAB} ${ML_CMD}

${fig_exp}: ${m_exp}
	${MATLAB} ${ML_CMD}

read:
	okular ${filename}.pdf &


clean:
	cd ${LATEX_ROOT} && make clean

clean_inkscape:
	cd ${LATEX_ROOT} && make clean_inkscape

clean_figures:
	cd ${LATEX_ROOT} && make clean_figures

clean_all:clean clean_inkscape clean_figures
