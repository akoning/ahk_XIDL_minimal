;WRAPPER FOR BUILT-IN HISTOGRAM FUNCTION THAT RETURNS A HISTOGRAM PLUS
;AN ARRAY OF X-AXIS VALUES FOR PLOTTING PURPOSES
function jdshist,array,min=min,max=max,binsize=binsize,outbins,plot=plot

;RUN HISTOGRAM
hist = histogram(array,min=min,max=max,binsize=binsize)

;MAKE X-AXIS
n = N_elements(hist)
outbins = min + binsize*findgen(n)

;PLOT (OR NOT)
if keyword_set(plot) then plot,outbins,hist,ps=10

return,hist
end
