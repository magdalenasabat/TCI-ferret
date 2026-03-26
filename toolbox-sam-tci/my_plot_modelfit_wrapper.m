function my_plot_modelfit_wrapper(M,L,folder_name)
% convinience wrapper for Sam's plotting function
% magdalena 14/02/2024s
    intper_sec = M.intper_sec;
    delay_sec_start = M.delay_sec_start;
    unique_segs = L.unique_segs;
    lag_t =  L.lag_t;
    intervaltype = 'highest';
    intervalmass = 0.75;

    plot_win = L.lag_t([1 end]); 
    plot_smoothwin = 0;
    linewidth = 2;

    plot_delaystat = 'median';
    plot_delay_range = [0, L.unique_segs(end)/1000];
    ploterrquant = 0.3;

    for chan=M.channels
        figh = figure;

        diff_context = L.diff_context(:,:,chan);
        same_context = L.same_context(:,:,chan);
        diff_context_bestpred  = M.diff_context_bestpred(:,:,chan);
        loss_best_model = M.loss(:,:,M.best_shape(chan)==M.shape, M.best_boundstrength(chan)==M.boundstrength, chan);
        best_intper_sec = M.best_intper_sec(chan);
        best_delay_sec_median = M.best_delay_sec_median(chan);
        best_shape = M.best_shape(chan);
        fname_global = mkpdir([folder_name '/ch-' num2str(chan) ]);


        plot_modelfit(diff_context, same_context, diff_context_bestpred, ...
            loss_best_model, best_intper_sec, best_delay_sec_median, best_shape, fname_global, ...
            intper_sec, delay_sec_start, unique_segs, lag_t,  intervaltype, intervalmass, ...
            plot_win, plot_smoothwin, plot_delaystat, plot_delay_range, ploterrquant, linewidth, figh)
    end
    close all
end