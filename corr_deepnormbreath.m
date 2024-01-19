clear;
close all;
dataDir='/media/yuhui/xxx';
cd(dataDir);


% group
normbh_list = [dir('xxx/swi_positive.sft/columnProfile/allcol1000.Breath.mean.vessels_outcount.bold_mdant.1D')];

for index=1:length(normbh_list)
	cd(dataDir);
	fprintf('++ Begin analyzing %s \n',normbh_list(index).name);

	subject_list = [dir(['*/swi_positive.sft/columnProfile/' normbh_list(index).name])];


	subjNum=length(subject_list);
	vesl0 = load([subject_list(1).folder '/' subject_list(1).name]);
	[~,columnNum] = size(vesl0);

	sub_normbreath_all = zeros(subjNum,columnNum);
	sub_deepbreath_all = zeros(subjNum,columnNum);


	for subj=1:subjNum
		cd(subject_list(subj).folder);

		effect0 = load(subject_list(subj).name);
		effect_norm=(effect0-0)/(mean(effect0(:))-0);
		sub_normbreath_all(subj,:) = effect_norm;
		fprintf('++ mean vesl = %f \n',mean(effect0(:)));

		if isfile('allcol1000.Breath.mean.act_rest.meanresponse.bold_mdant.1D')

			effect0 = load('allcol1000.Breath.mean.act_rest.meanresponse.bold_mdant.1D');
			effect_norm = (effect0-0)/(mean(effect0(:))-0);
			sub_deepbreath_all(subj,:) = effect_norm;
			fprintf('++ mean depbreath = %f \n',mean(effect0(:)));
		end

	end

	cd([dataDir '/group_allcol1000']);

	%% sub_deepbreath_all .................................................
	temps = mean(sub_deepbreath_all,2);
	goodsubj = find(temps~=0);

	sub_deepbreath_all_gd = sub_deepbreath_all(goodsubj,:);
	sub_normbreath_all_gd = sub_normbreath_all(goodsubj,:);

	[r,p] = corrcoef(sub_normbreath_all_gd(:),sub_deepbreath_all_gd(:));

	mdl = fitlm(sub_normbreath_all_gd(:),sub_deepbreath_all_gd(:));


	xfit=sub_normbreath_all_gd(:);
	yfit=sub_deepbreath_all_gd(:);
	figure;
	h=plot(mdl);
	% Get handles to plot components
	dataHandle = findobj(h,'DisplayName','Data');
	fitHandle = findobj(h,'DisplayName','Fit');
	% The confidence bounds have 2 handles but only one of 
	% the handles contains the legend string.  The first
	% line below finds that object and then searches for 
	% other objects in the plot that have the same linestyle
	% and color. 
	cbHandles = findobj(h,'DisplayName','Confidence bounds');
	cbHandles = findobj(h,'LineStyle',cbHandles.LineStyle, 'Color', cbHandles.Color);

	dataHandle.Color = 'none'; 
	% fitHandle.Color = [1 0 0]; %orange 
	set(fitHandle, 'Color', 'red', 'LineWidth', 1.2)
	set(cbHandles, 'Color', 'red', 'LineWidth', 1.0)
	legend('off');

	ylabel('Deep breath','Fontsize',20,'FontWeight','bold');
	xlabel('Natural breath','Fontsize',20,'FontWeight','bold');

	hold on;
	histogram2(xfit,yfit,100,'DisplayStyle','tile','EdgeColor','none','ShowEmptyBins','off','FaceAlpha',1); % 
	colormap(flipud(winter)); % bone gray pink
	grid off
	% colorbar;
	set(gca,'CLim',[3 25]);
	box off
	xlim([0 3]);
	ylim([0 3]);
	whitebg('white');
	set(gcf,'color',[1 1 1]);
	set(gca,'linewidth',3,'fontsize',20,'FontWeight','normal','Xcolor',[0 0 0],'Ycolor',[0 0 0])
	% title(['r = ' num2str(r(1,2)) ', p = ' num2str(p(1,2))],'fontsize',18,'FontWeight','normal');
	title(['r = ' num2str(r(1,2))],'fontsize',18,'FontWeight','normal');
	export_fig([extractBefore(normbh_list(index).name,'.1D') '.sub_nomdepbreath.png'],'-r300');
end





