clear all
clc


set(gcf,'position',[0,0,700,700]);


datasets = {'DHFR', 'PROTEINS_full'};


for idata = 1:length(datasets)
    dataset = datasets{idata};
    load(['results/',dataset, '/', dataset, '_hashwls_results.mat']);

    subplot(2,2, idata)
    
    plot(1:5, accs_mean(:,1), '-+','Color',[0.9412 0.4706 0],'LineWidth',1)
    hold on;
    plot(1:5, accs_mean(:,2), '-rx','LineWidth',1)
    plot(1:5, accs_mean(:,3), '-d','Color',[0.502 0.251 0],'LineWidth',1)
    plot(1:5, accs_mean(:,4), '-*','Color',[0 0.502 0.502],'LineWidth',1)
    plot(1:5, accs_mean(:,5), '-o','Color',[0.502 0.502 0.502],'LineWidth',1)
    plot(1:5, accs_mean(:,6), '-.^b','LineWidth',1)

    legend({'$K=50$', '$K=100$', '$K=150$', '$K=200$', '$K=250$', '$K=300$'}, 'location', 'northwest', 'Interpreter', 'latex', 'FontSize', 10)
    set(legend,'color','none')
    xlabel('Iterations', 'FontSize', 14)
    set(gca, 'xtick', [1 2 3 4 5], 'xticklabels', {'1', '2', '3', '4', '5'})
    ylabel('Accuracy (%)', 'FontSize', 14)
    xlim([0 6])
    if idata == 1
        ylim([65,85])
    elseif idata == 2
        ylim([65,85])
    end
    title([dataset], 'FontSize', 14,'Interpreter','none')
    
    subplot(2,2, idata+2)
    
    plot(1:5, cpus_mean(:,1), '-+','Color',[0.9412 0.4706 0],'LineWidth',1)
    hold on;
    plot(1:5, cpus_mean(:,2), '-rx','LineWidth',1)
    plot(1:5, cpus_mean(:,3), '-d','Color',[0.502 0.251 0],'LineWidth',1)
    plot(1:5, cpus_mean(:,4), '-*','Color',[0 0.502 0.502],'LineWidth',1)
    plot(1:5, cpus_mean(:,5), '-o','Color',[0.502 0.502 0.502],'LineWidth',1)
    plot(1:5, cpus_mean(:,6), '-.^b','LineWidth',1)
    legend({'$K=50$', '$K=100$', '$K=150$', '$K=200$', '$K=250$', '$K=300$'}, 'location', 'northwest', 'Interpreter', 'latex', 'FontSize', 10)
    set(legend,'color','none')
    xlabel('Iterations', 'FontSize', 14)
    set(gca, 'xtick', [1 2 3 4 5], 'xticklabels', {'1', '2', '3', '4', '5'})
    ylabel('Runtime (s)', 'FontSize', 14)
    xlim([0 6])
    if idata == 1
        ylim([0,30])
    elseif idata==2
        ylim([0,50])
    end
    title([dataset], 'FontSize', 14,'Interpreter','none')
    
end
