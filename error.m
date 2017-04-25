% function [outcoef,se] = error(coef,data,indvars,jacobain)
% 
% 
% 	startfit = 1; 
%     endfit = 100;
% 	fit = [];
%     length=endfit-startfit;
% 	nindvars = length(indvars(:,1));
% 	npts = length(indvars);  % number of channels
%     invm = []; 
%     outcoef= [];
%     invm = inv(jacobain);     
% 			for i=1:1:npts
% 				for j=1:1:nindvar
% 					fit(i) = fit(i) + coef(j)*indvars(j,i);
%                 end
%             end
%             
%             	
% 			for i = startfit:1:endfit
% 				tempc2= tempc2 + (fit(i)-data(i))*(fit(i)-data(i));
%             end
%             
%             c2 = tempc2/(length-nindvars);
% 
% 
%                 se = sqrt(c2*invm(i,i)); %std?
% 
% 		for i=1:1:length(invm)
% 			for j=1:1:length(jvector)
% 				outcoef(i)= outcoef(i) + invm(i,j)*jvector(j);
%             end
%         end
% 		
% end
% 
% 
%             
            
            
% 		double[] se=new double[nindvars];
% 		double c2=get_c2(outcoef,data,weights);
% 		//int length=endfit-startfit+1;
% 		//c2*=(double)(length-nindvars)/(double)(length-1);
% 		for(int i=0;i<nindvars;i++) se[i]=Math.sqrt(c2*inv[i][i]);
% 		return new double[][]{outcoef,se};
%         
%         
%         public double get_c2(double[] coef,float[] data,float[] weights){
% 		//error

