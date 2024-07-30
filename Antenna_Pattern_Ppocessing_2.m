classdef Antenna_Pattern_Ppocessing_2 < matlab.apps.AppBase

    % Properties that correspond to app components
    properties (Access = public)
        UIFigure                  matlab.ui.Figure
        Grid                      matlab.ui.container.GridLayout
        DropDown_step             matlab.ui.control.DropDown
        Format                    matlab.ui.control.DropDown
        FormatLabel               matlab.ui.control.Label
        IncidentWaveARRwPLFLabel  matlab.ui.control.Label
        LossindBLabel             matlab.ui.control.Label
        Button_Coverage           matlab.ui.control.Button
        DropDown                  matlab.ui.control.DropDown
        Switch                    matlab.ui.control.Switch
        Export_Output             matlab.ui.control.Button
        Export_Input              matlab.ui.control.Button
        Table_DataOut             matlab.ui.control.Table
        Table_DataIn              matlab.ui.control.Table
        Label_Pol                 matlab.ui.control.Label
        Button_ContourPlot        matlab.ui.control.Button
        Button_3dPlot             matlab.ui.control.Button
        Label_max                 matlab.ui.control.Label
        Panel_Rect                matlab.ui.container.Panel
        Grid_Rect                 matlab.ui.container.GridLayout
        Axes2                     matlab.ui.control.UIAxes
        EditField_Rw              matlab.ui.control.Spinner
        Panel_Polar               matlab.ui.container.Panel
        Grid_Polar                matlab.ui.container.GridLayout
        GridLayout                matlab.ui.container.GridLayout
        CheckBox_Et               matlab.ui.control.CheckBox
        CheckBox_Er               matlab.ui.control.CheckBox
        CheckBox_El               matlab.ui.control.CheckBox
        Button_ExportCut          matlab.ui.control.Button
        RangeMax                  matlab.ui.control.Spinner
        RangeMin                  matlab.ui.control.Spinner
        Label_HPBW                matlab.ui.control.Label
        Button_HPBW               matlab.ui.control.StateButton
        Range                     matlab.ui.control.RangeSlider
        EditField_Loss            matlab.ui.control.Spinner
        Button_Load               matlab.ui.control.Button
    end

    
    properties (Access = private)
        FilePath char % InputFile Name
        FolderPath char % InputFile Name
        InputFileName char
        Rmax double
        Rmin double
        Phase double
        Gmax double
        Gmax_Theta double
        step double
        DataIn_step table
        DataOut_step table
        cut_step table
        cutE table
        cutH table
        cut table
        Theta_Phi_0_180 table
        Theta_Phi_90_27 table
        Phi_Theta90 table
        pax matlab.graphics.axis.PolarAxes
        Plot1 matlab.graphics.chart.primitive.Line
        Plot1_2 matlab.graphics.chart.primitive.Line
        Plot1_3 matlab.graphics.chart.primitive.Line
        maxGain_RLine matlab.graphics.chart.primitive.Line
        maxGain_XLine matlab.graphics.primitive.Line
        maxGain_YLine matlab.graphics.primitive.Line
        region_Polar matlab.graphics.chart.decoration.PolarRegion
        region_Rect  matlab.graphics.chart.decoration.ConstantRegion
        % ButtonHPWB matlab.ui.control.StateButton
        % Button matlab.ui.control.StateButton
    end
    
    methods (Access = private)
        
        function loadData(app, ~)
            app.Switch.Visible = 'on';
            app.EditField_Loss.Visible = 'on';
            app.LossindBLabel.Visible = 'on';
            app.Export_Input.Visible = 'on';
            app.Export_Output.Visible = 'on';
            app.Button_3dPlot.Visible = 'on';
            app.Button_ContourPlot.Visible = 'on';
            app.Button_Coverage.Visible = 'on';
            app.Table_DataIn.Visible = 'on';
            app.Table_DataOut.Visible = 'on';
            app.Panel_Polar.Visible = 'on';
            app.Panel_Rect.Visible = 'on';
            app.DropDown.Visible = 'on';
            app.EditField_Rw.Visible  = 'on';
            app.IncidentWaveARRwPLFLabel.Visible  = 'on';
            app.CheckBox_Et.Value = true;
            app.CheckBox_Er.Value = true;
            app.CheckBox_El.Value = true;

           % Detect headerLines
            headerLines = max(0, find(~cellfun(@isempty, regexp(strtrim(readlines(app.FilePath)), '^[+-]?\d+(\.\d+)?([eE][+-]?\d+)?', 'once')), 1) - 1);
            
            % Read the data
            opts = detectImportOptions(app.FilePath, 'FileType','text','NumHeaderLines',headerLines,'LeadingDelimitersRule','ignore','ConsecutiveDelimitersRule','join');
            opts.VariableNames = {'Theta', 'Phi', 'E_TH_DB', 'E_PH_DB', 'E_TH_DG', 'E_PH_DG'};
            dataIn = readtable(app.FilePath,opts);
            
            % Set loss/offset
            % if contains(app.FilePath, {'PPE','nX','-X'}, 'IgnoreCase', true)
            if contains(app.FilePath, 'PPE', 'IgnoreCase', true)
                app.EditField_Loss.Value = 0;
            end
            % loss{i} = -Lc - (Lp - Lc) * contains(subfolders2{i}, {'PPE','nX','-X'}, 'IgnoreCase', true);
            loss = app.EditField_Loss.Value;
            dataIn.E_TH_DB = dataIn.E_TH_DB - loss;
            dataIn.E_PH_DB = dataIn.E_PH_DB - loss;

            app.Table_DataIn.Data = dataIn;
            app.DataIn_step = dataIn;

            Theta_unique = unique(dataIn.Theta);
            app.step = Theta_unique(2) - Theta_unique(1);
            app.DropDown_step.Items = {['STEP: ', num2str(app.step), '°'], 'STEP: 1°'};
            if app.step < 1
                app.DropDown_step.Enable = 'on';
                app.DropDown_step.Visible = 'on';
                app.DropDown_step.Placeholder = 'STEP';
                % app.DropDown_stepValueChanged;
            else
                app.DropDown_step.Value = ['STEP: ', num2str(app.step), '°'];
                app.DropDown_step.Visible = 'off';
                app.DropDown_step.Enable = 'off';
            end

            app.processData();
        end

        function processData(app)
            if app.step < 1 && isequal(app.DropDown_step.Value,'STEP: 1°')
                app.Table_DataIn.Data = app.Table_DataIn.Data(ismember(app.Table_DataIn.Data.Theta, (0:180)) & ismember(app.Table_DataIn.Data.Phi, (0:360)),:);
                % app.Table_DataOut.Data = app.Table_DataOut.Data(ismember(app.Table_DataOut.Data.Theta, (0:180)) & ismember(app.Table_DataOut.Data.Phi, (0:360)),:);
                % app.cut = app.cut(ismember(app.cut.Theta, (0:360)) & ismember(app.cut.Phi, (0:360)),:);
            else
                app.Table_DataIn.Data = app.DataIn_step;
                % app.Table_DataOut.Data = app.DataOut_step;
                % app.cut = app.cut_step;
            end

            dataIn = app.Table_DataIn.Data;
            if app.Format.Value == 1 % Format1: Theta Phi (mag_dB, phase°)
                dataIn.Properties.VariableNames = {'Theta', 'Phi', 'E_TH_DB', 'E_PH_DB', 'E_TH_DG', 'E_PH_DG'};
                app.Table_DataIn.ColumnName = {'Theta', 'Phi', 'E_TH_DB', 'E_TH_DB', 'E_TH_DG', 'E_TH_DG'};
                POL1_dB = dataIn.E_TH_DB;
                POL2_dB = dataIn.E_PH_DB;

                % Compute E_theta & E_phi in complex form
                Eth = complex(db2mag(dataIn.E_TH_DB) .* cosd(dataIn.E_TH_DG) , db2mag(dataIn.E_TH_DB) .* sind(dataIn.E_TH_DG));
                Eph = complex(db2mag(dataIn.E_PH_DB) .* cosd(dataIn.E_PH_DG) , db2mag(dataIn.E_PH_DB) .* sind(dataIn.E_PH_DG));
                
                % Compute E_rcp, E_lcp
                Ercp = (Eth + 1i*Eph)./sqrt(2);
                Elcp = (Eth - 1i*Eph)./sqrt(2);
                Ercp_mag = abs(Ercp);
                Elcp_mag = abs(Elcp);
            elseif app.Format.Value == 2 % Format2: Theta Phi (Re, Im)
                dataIn.Properties.VariableNames = {'Theta', 'Phi', 'E_TH_Re', 'E_TH_Im', 'E_PH_Re', 'E_PH_Im'};
                app.Table_DataIn.ColumnName = {'Theta', 'Phi', 'E_TH_Re', 'E_TH_Im', 'E_PH_Re', 'E_PH_Im'};

                POL1_dB = 10.^(dataIn.E_TH_Re/20);
                POL2_dB = 10.^(dataIn.E_PH_Re/20);

                % Compute E_theta & E_phi in complex form
                Eth = complex(dataIn.E_TH_Re , dataIn.E_TH_Im);
                Eph = complex(dataIn.E_PH_Re , dataIn.E_PH_Im);
               
                % Compute E_rcp, E_lcp
                Ercp = (Eth + 1i*Eph)./sqrt(2);
                Elcp = (Eth - 1i*Eph)./sqrt(2);
                Ercp_mag = abs(Ercp);
                Elcp_mag = abs(Elcp);
            elseif app.Format.Value == 3 % Format3: Ercp Eclp (mag, phase) - POL1=RCP POL2=LCP
                dataIn.Properties.VariableNames = {'Theta', 'Phi', 'POL1_DB', 'POL2_DB', 'POL1_DG', 'POL2_DG'};
                app.Table_DataIn.ColumnName = {'Theta', 'Phi', 'E_TH_DB', 'E_TH_DB', 'E_TH_DG', 'E_TH_DG'};

                Ercp_mag = 10.^(POL1_DB/20);
                Elcp_mag = 10.^(POL2_DB/20);
                % Ercp_mag = 20.*log10(POL1_DB);
                % Elcp_mag = 20.*log10(POL2_DB);

            else % Format4: Ercp Eclp (Re, Im)
                dataIn.Properties.VariableNames = {'Theta', 'Phi', 'POL1_real', 'POL1_imag', 'POL2_real', 'POL2_imag'};
                app.Table_DataIn.ColumnName = {'Theta', 'Phi', 'POL1_real', 'POL1_imag', 'POL2_real', 'POL2_imag'};
                % set(Table, 'ColumnName', {'Theta', 'Phi', 'Pol1', 'Pol2'});
                % POL1_dB = dataIn.POL1_real;
                % POL2_dB = dataIn.POL2_real;

                % Compute E_rcp, E_lcp
                Ercp_mag = sqrt(dataIn.POL1_real.^2 + dataIn.POL1_imag.^2);
                Elcp_mag = sqrt(dataIn.POL2_real.^2 + dataIn.POL2_imag.^2);

                POL1_dB = 20.*log10(Ercp_mag);
                POL2_dB = 20.*log10(Elcp_mag);
            end


            % Compute Total Gain
            E_total = sqrt(Ercp_mag.^2+Elcp_mag.^2);
            
            % Compute Axial Ratio
            AR = (Ercp_mag + Elcp_mag)./(Ercp_mag - Elcp_mag);
            AR (Ercp_mag == Elcp_mag) = sqrt(10)/1e+13;

            sgn = sign(Ercp_mag - Elcp_mag);
            sgn (Ercp_mag == Elcp_mag) = -1;
            % AR = sgn.*(20*log10( (Ercp_mag + Elcp_mag)./abs(Ercp_mag - Elcp_mag) ));
            % AR (Ercp_mag == Elcp_mag) = -999; % AR (AR == Inf) =  -999;

            % Compute Polarization Loss
            tau = repmat(90, size(AR)); % relative titl angle
            % tau (AR<0) = 0;
            Ra = AR; % 10.^(abs(AR)./20); % axial ratio of Antenna (Receiving Antenna 'a')    %dB2mag
            % Rw = 10.^(repmat(app.EditField_Rw.Value, size(Ra))./20); % axial ratio of Wave (Incident Wave 'w')
            Rw_db = app.EditField_Rw.Value;
            Rw = repmat(10.^(Rw_db/20), size(AR)); % axial ratio of Wave (Incident Wave 'w')
            % PLF = 10*log10(0.5 +  sgn.*(4*Ra.*Rw + (Ra.^2 - 1).*(Rw.^2 - 1).*cosd(2*tau))./(2*(Ra.^2 + 1).*(Rw.^2 + 1)));
            PLF = 10*log10(0.5 +  (4*Ra.*Rw + (Ra.^2 - 1).*(Rw.^2 - 1).*cosd(2*tau))./(2*(Ra.^2 + 1).*(Rw.^2 + 1)));
            % additonal optional Polarization Loss cases with Rw = 2dB & 7dB

        
            % Gain data table (Total Gain, CP components & corrected gain)
            dataOut = dataIn;
            varNames = {'Theta', 'Phi', 'E_Total_dB', 'E_RCP_dB', 'E_LCP_dB', 'AR_dB'};
            dataOut.Properties.VariableNames = varNames;
            dataOut.E_Total_dB = 20.*log10(E_total);
            dataOut.E_RCP_dB = 20.*log10(Ercp_mag);        % RCP_mag_db=20.*log10(Ercp_mag) | RCP_Phase = rad2deg(angle(Ercp));
            dataOut.E_LCP_dB = 20.*log10(Elcp_mag);        % LCP_mag_db=20.*log10(Elcp_mag) | LCP_Phase = rad2deg(angle(Elcp));
            dataOut.AR_dB = 20.*log10(abs(AR)).*sign(AR);  % Axial Ratio=20.*log10(mag)
            dataOut.PLF = PLF;                          % Polarization Loss
            dataOut.Pol_Gain = dataOut.E_Total_dB + PLF; % Generic VV AR=6dB

            app.Table_DataOut.Data = dataOut;
            app.DataOut_step = dataOut;

            % Display Max Theta-Phi Component
            % if gt(max(dataIn.E_TH_DB),max(dataIn.E_PH_DB))
            %     mx=max(dataIn.E_TH_DB); component=' dB (E-TH)';
            % elseif lt(max(dataIn.E_TH_DB),max(dataIn.E_PH_DB))
            %     mx=max(dataIn.E_PH_DB); component=' dB (E-PH)';
            % else
            %     mx=max(dataIn.E_TH_DB); component=' dB';
            % end
            
            if gt(max(POL1_dB),max(POL2_dB))
                mx=max(POL1_dB); component=' dB (E-TH)';
            elseif lt(max(POL1_dB),max(POL2_dB))
                mx=max(POL2_dB); component=' dB (E-PH)';
            else
                mx=max(POL1_dB); component=' dB';
            end

            app.Label_max.Visible = 'on';
            app.Label_max.Text = ['Max E_{θ/φ} component (xGTD Max-Gain Input):  ' num2str(mx) component];
            % app.Label_max.Text = ['Maximum Theta/Phi component: ' num2str(mx) component];

            POB = max(dataOut.E_Total_dB);

            % Display Antenna Pattern Polarization
            if gt(sum(dataOut.E_RCP_dB),sum(dataOut.E_LCP_dB)), pol='RCP-Polarized'; else, pol='LCP-Polarized'; end
            % delta = sum(dataOut.E_RCP_dB) - sum(dataOut.E_LCP_dB);

            app.Label_Pol.Visible = 'on';
            % app.Label_Pol.Text = pol;
            app.Label_Pol.Text = [pol , ' | POB:  ' num2str(round(POB,2)) 'dBi'];



            assignin('base', 'step', app.step);

            % Get E-Plane & H-Plane cuts
            % Theta_Phi_0_180 cut
            Theta_Phi0 = app.Table_DataOut.Data(app.Table_DataOut.Data.Phi ==  0, :);
            Theta_Phi180 = flip(app.Table_DataOut.Data(app.Table_DataOut.Data.Phi ==  180, :));
            Theta_Phi180.Theta = flip(Theta_Phi180.Theta) + 180; % = 360 - Theta_Phi180.Theta;
            app.Theta_Phi_0_180 = vertcat(Theta_Phi0 , Theta_Phi180(2:end,:));
            % Theta_Phi_90_270 cut
            Theta_Phi90 = app.Table_DataOut.Data(app.Table_DataOut.Data.Phi ==  90, :);
            Theta_Phi270 = flip(app.Table_DataOut.Data(app.Table_DataOut.Data.Phi ==  270, :));
            Theta_Phi270.Theta = flip(Theta_Phi270.Theta) + 180;
            app.Theta_Phi_90_27 = vertcat(Theta_Phi90 , Theta_Phi270(2:end,:));
            % Phi_Theta90 cut
            app.Phi_Theta90 = app.Table_DataOut.Data(app.Table_DataOut.Data.Theta == 90, :);
            
            switch true
                case contains(app.FilePath, {'pZ','Zp','+Z','Z+'}, 'IgnoreCase', true),  orientation = '+Z';
                case contains(app.FilePath, {'nZ','Zn','-Z','Z-'}, 'IgnoreCase', true),  orientation = '-Z';
                case contains(app.FilePath, {'pX','Xp','+X','X+'}, 'IgnoreCase', true),  orientation = '+X';
                case contains(app.FilePath, {'nX','Xn','-X','X-'}, 'IgnoreCase', true),  orientation = '-X';
                case contains(app.FilePath, {'pY','Yp','+Y','Y+'}, 'IgnoreCase', true),  orientation = '+Y';
                case contains(app.FilePath, {'nY','Yn','-Y','Y-'}, 'IgnoreCase', true),  orientation = '-Y';
                otherwise
                    orientation = 'Orientation:';
            end
            app.DropDown.Value = orientation;
            app.DropDownValueChanged();
        end
        
        function plotPattern(app, ~)

            % if app.Switch.Value == 'E-Plane' (% isequal ==)
            app.cut = app.cutE;
            app.Phase = app.cut.Theta;
            cutAngle = 'Theta (degree)';
            
            if app.Switch.Value == "H-Plane"
                app.cut = app.cutH;
                % if ismember(app.DropDown.Value, ["X", "Y"]) %if app.DropDown.Value == "X-axis" || app.DropDown.Value == "Y-axis"
                % if ismember({app.DropDown.Value}, {'X', 'Y'}) %if app.DropDown.Value == "X-axis" || app.DropDown.Value == "Y-axis"
                if contains(app.DropDown.Value, {'X', 'Y'})
                    app.Phase = app.cut.Phi;
                    cutAngle = 'Phi (degree)';
                end
            end

            app.cut_step = app.cut;

            [app.Gmax, maxGain_idx] = max(app.cut.E_Total_dB);
            app.Gmax_Theta = app.Phase(maxGain_idx);

            xlabel(app.Axes2, cutAngle);

            % pattern = app.cut;
            Gtot = app.cut.E_Total_dB;
            Grcp = app.cut.E_RCP_dB;
            Glcp = app.cut.E_LCP_dB;

            epsilon = 0.5; % to adjust the rounding
            app.Rmin = floor(min(min([Gtot,Grcp,Glcp]))/10 - epsilon)*10;
            if app.Rmin < -50
                app.Rmin = -50;
            end
            
            app.Rmax = ceil(max(max([Gtot,Grcp,Glcp]))/10 + epsilon)*10;

            app.RangeMin.Value = app.Rmin;
            app.RangeMax.Value = app.Rmax;

            app.Range.Limits = [app.Rmin app.Rmax];
            app.Range.Value = [app.Rmin app.Rmax];
            MajorTicksStep = 10;
            app.Range.MajorTicks = app.Rmin:MajorTicksStep:app.Rmax;
            app.Range.MajorTickLabels = app.Range.MajorTicks + " dB";

            % if ishandle(app.pax.Children)
            if isgraphics(app.pax.Children)
                delete(app.pax.Children);
            end

            app.Plot1 = polarplot(app.pax,app.Phase*pi/180,Gtot,'LineWidth',2,'Color','red',DisplayName='E_Total');
            % title(pax,'Polar Plot');
            thetaticks(app.pax,0:15:360);

            app.pax.RLim = [app.Rmin app.Rmax]; %rlim(app.pax,[app.Rmin app.Rmax]);
            app.pax.RTick = (app.Rmin:10:app.Rmax); %rticks(app.pax,(app.Rmin:10:app.Rmax))
            app.pax.ThetaZeroLocation='top';
            app.pax.ThetaDir = 'clockwise';
            hold(app.pax, 'on');
            app.Plot1_2 = polarplot(app.pax,app.Phase*pi/180,Grcp,'LineWidth',2,'Color','#0072BD',DisplayName='E_RCP');
            app.Plot1_3 = polarplot(app.pax,app.Phase*pi/180,Glcp,'LineWidth',2,'Color','#77AC30',DisplayName='E_LCP');
            % legend('show');
            legend(app.pax,[app.Plot1,app.Plot1_2,app.Plot1_3], 'Interpreter', 'none', 'Location','northeastoutside');
            % uistack(Plot1_2);
            uistack(app.Plot1,'top');
            hold(app.pax, 'off');

            % Rect Plot
            % cla(app.Axes2);
            app.Axes2.Visible = 'on';
            Plot2 = plot(app.Axes2,app.Phase,Gtot, 'LineWidth', 2, 'Color', 'red',DisplayName='E_Total');
            % a = annotation('textarrow',app.Gmax,app.Gmax_Theta,'String',['Gain_max  = ' num2str(round(app.Gmax,2)) 'dB']);
            % idxmax = find(Gtot == max(Gtot));
            % Plot2 = plot(app.Axes2,app.Phase,Gtot, '-o' ,'MarkerIndices', idxmax, 'MarkerEdgeColor','k', 'MarkerFaceColor','k', 'MarkerSize',2.5, 'LineWidth', 2, 'Color', 'red', DisplayName='E_Total');
            ylim(app.Axes2,[app.Rmin app.Rmax]);
            yticks(app.Axes2,(app.Rmin:5:app.Rmax));
            % xlim(app.Axes2,[0 180]);
            grid(app.Axes2,'on');
            % x=app.Phase(1);
            % y=max(Gtot);
            x=app.Gmax_Theta;
            y=app.Gmax;
            % disp(['X = ' , num2str(x)] );
            % disp(['y = ' , num2str(y)] );
            % disp(y);
            % x=pattern.Theta(1);
            % y=pattern.E_Total_dB(1);
            % annotation(app.Axes2,'textarrow',x,y,'String',['Boregiht: ' num2str(y) 'dB']);
            % text(app.Axes2,x,max(app.Axes2.YTick),['Gain_{MAX} = ' num2str(round(y,2)) 'dB'], 'Units','normalized', 'Position', [0.01 0.95]);
            % text(app.Axes2,app.Gmax_Theta,app.Gmax,['Gain_{MAX} = ' num2str(round(app.Gmax,2)) 'dB']);
            % text(app.Axes2,app.Phase(1),max(app.Axes2.YTick),['Gain_{MAX} = ' num2str(round(y,2)) 'dB'], 'Units','normalized', 'Position', [0.01 0.95]);
            max_text = text(app.Axes2,x,y,['Gain_{MAX} = ' num2str(round(y,2)) 'dB'],'VerticalAlignment','bottom','FontWeight','bold');
            X_line = line(app.Axes2, [x x], [app.Rmin y], 'LineStyle', '--', 'Color', 'k');
            Y_line = line(app.Axes2, [app.Rmin x], [y y], 'LineStyle', '--', 'Color', 'k');
            % text(app.Axes2,app.Gmax,max(app.Axes2.YTick),['Gain_{MAX} = ' num2str(round(app.Gmax,2)) 'dB'], 'Units','normalized', 'Position', [0.01 0.95]);
            hold(app.Axes2, 'on');
            Plot2_2 = plot(app.Axes2,app.Phase,Grcp, 'LineWidth', 1.5, 'Color', '#0072BD',DisplayName='E_RCP');
            Plot2_3 = plot(app.Axes2,app.Phase,Glcp, 'LineWidth', 1.5, 'Color', '#77AC30',DisplayName='E_LCP');
            % xlim(app.Axes2,[0 180]); grid(app.Axes2,'on');
            legend(app.Axes2,[Plot2,Plot2_2,Plot2_3], 'Interpreter', 'none', 'Location','northeastoutside');
            uistack(Plot2,'top');
            hold(app.Axes2, 'off');

            if app.Button_HPBW.Value
                app.Button_HPBWValueChanged();
            end
        end
    end
    

    % Callbacks that handle component events
    methods (Access = private)

        % Code that executes after component creation
        function startupFcn(app)
            app.Table_DataIn.RowName = 'numbered';
            app.Table_DataOut.RowName = 'numbered';

            app.pax = polaraxes(app.Grid_Polar);
            cla(app.pax,"reset");
            reset(app.pax);
            app.pax.Layout.Row = [1 4];
            app.pax.Layout.Column = 3;

            app.cut = app.cutE;
        end

        % Button pushed function: Button_Load
        function Button_LoadPushed(app, event)
            [FileName, app.FolderPath] = uigetfile({'*.uan;*.fz;*.cut;*.out;*.csv;*.txt','Antenna Pattern file'}, 'Select an Antenna Pattern File');
            if isequal(FileName,0) || isequal(app.FolderPath,0)
                return; % User canceled the operation
            else
                app.FilePath = fullfile(app.FolderPath, FileName);
                [~, app.InputFileName, ~] = fileparts(app.FilePath);
                app.UIFigure.Name = ['PatternReader - ' app.FilePath];
                app.loadData();
            end
        end

        % Value changed function: EditField_Loss
        function EditField_LossValueChanged2(app, event)
            if ~isempty(app.Table_DataIn.Data)
                app.processData();
            end
        end

        % Button pushed function: Export_Input
        function Export_InputButtonPushed(app, event)
            % default = regexp(app.FilePath;
            if ~isempty(app.Table_DataIn.Data)
                cd(app.FolderPath);
                [FileName, PathName, filterindex] = uiputfile({'*.csv', 'CSV Files (*.csv)'; '*.txt', 'Text Files (*.txt)'}, 'Save Input Data',[app.FilePath, '_input_', char(extractBetween(app.DropDown_step.Value, 'STEP: ', '°')), '_deg']);
                if isequal(FileName,0) || isequal(PathName,0)
                    return; % User canceled the operation
                else
                    filepath = fullfile(PathName, FileName);
                    if filterindex == 1
                        % Export as CSV
                        writetable(app.Table_DataIn.Data, filepath);
                    elseif filterindex == 2
                        % Export as TXT
                        writetable(app.Table_DataIn.Data, filepath, 'Delimiter', '\t');
                    end
                end
            else
                msgbox('No data to export.', 'Warning', 'warn');
            end
        end

        % Button pushed function: Export_Output
        function Export_OutputButtonPushed(app, event)
            if ~isempty(app.Table_DataOut.Data)
                cd(app.FolderPath);
                [FileName, PathName, filterindex] = uiputfile({'*.csv', 'CSV Files (*.csv)'; '*.txt', 'Text Files (*.txt)'}, 'Save Input Data',[app.InputFileName, '_output_', char(extractBetween(app.DropDown_step.Value, 'STEP: ', '°')) , '_deg']);
                if isequal(FileName,0) || isequal(PathName,0)
                    return; % User canceled the operation
                else
                    filepath = fullfile(PathName, FileName);
                    if filterindex == 1
                        % Export as CSV
                        writetable(app.Table_DataOut.Data, filepath);
                    elseif filterindex == 2
                        % Export as TXT
                        writetable(app.Table_DataOut.Data, filepath, 'Delimiter', '\t');
                    end
                end
            else
                msgbox('No data to export.', 'Warning', 'warn');
            end
        end

        % Button pushed function: Button_3dPlot
        function Button_3dPlotPushed(app, event)
            % Data = app.Table_DataOut.Data;
            % fig_3D = gcf;
            % fig_3d.Name = '3D Plot';
            figure('Name','3D Plot','NumberTitle','off');
            plot3D = patternCustom(app.Table_DataOut.Data.E_Total_dB,app.Table_DataOut.Data.Theta,app.Table_DataOut.Data.Phi); %3D
        end

        % Value changed function: EditField_Rw
        function EditField_RwValueChanged2(app, event)
            if ~isempty(app.Table_DataIn.Data)
                app.processData();
            end
        end

        % Value changed function: Range
        function RangeValueChanged(app, event)
            % dbScale = app.Range  .Value;
            app.Rmin = min(app.Range.Value);
            app.Rmax = max(app.Range.Value);
            rlim(app.pax,[app.Rmin app.Rmax]);
            ylim(app.Axes2,[app.Rmin app.Rmax]);
            % yticks(app.Axes2,(app.Rmin:5:app.Rmax));
        end

        % Value changed function: Button_HPBW
        function Button_HPBWValueChanged(app, event)
            % value = app.Button_HPBW.Value;
            if app.Button_HPBW.Value
                app.Button_HPBW.BackgroundColor = [0 1 0];  % Green

                gain = app.cut.E_Total_dB;

                [maxGain, maxGain_Idx] = max(gain);
                app.Gmax = maxGain;
                Gain3dB = maxGain - 3;

                theta_maxGain = app.Phase(maxGain_Idx);
                hold(app.pax, 'on');
                app.maxGain_RLine = polarplot(app.pax, [theta_maxGain theta_maxGain]*pi/180, [app.Rmin maxGain], 'k--', 'DisplayName', 'Max Gain');
                app.maxGain_RLine.Annotation.LegendInformation.IconDisplayStyle = "off";
                hold(app.pax, 'off');
                
                app.Gmax_Theta = theta_maxGain;

                % binnary array
                Gain_above3dB_bin = ones(size(gain));
                Gain_above3dB_bin(gain < Gain3dB) = 0;

                boundary_Right_Idx = find(Gain_above3dB_bin(maxGain_Idx:end) == 0,1,'first') - 2 + maxGain_Idx;
                boundary_Left_Idx = find(Gain_above3dB_bin(1:maxGain_Idx) == 0,1,'last') + 1;

                HP_array = Gain_above3dB_bin;

                if isempty(boundary_Left_Idx)
                    crossZero = true;
                    boundary_Left_Idx = find(Gain_above3dB_bin == 0,1,'last') + 1;
                    HP_array (boundary_Right_Idx:boundary_Left_Idx) = 0;
                else
                    crossZero = false;
                    HP_array (1:boundary_Left_Idx-1) = 0;
                    HP_array (boundary_Right_Idx+1:end) = 0;
                end
                
                Theta = app.Phase;

                
                HPBW_array = nonzeros(HP_array);
                HPBW = numel(HPBW_array)*app.step;
                
                theta_bounds = Theta([boundary_Right_Idx,boundary_Left_Idx]);
                
                if crossZero
                    theta_bounds = wrapTo180(theta_bounds);
                end
                
                % HPBW_region
                app.region_Polar = thetaregion(app.pax, theta_bounds*pi/180, FaceColor='y', FaceAlpha=0.3, DisplayName='HPBW');
                % app.region_Polar = thetaregion(app.pax, theta_bounds*pi/180, FaceColor='y', FaceAlpha=0.3, EdgeColor='#EDB120', LineStyle='--', DisplayName='HPBW');
                app.region_Polar.Annotation.LegendInformation.IconDisplayStyle = "off";
                app.region_Rect = xregion(app.Axes2 , theta_bounds, FaceColor='y', DisplayName='HPBW'); %FaceColor='#a9ddff', EdgeColor='#EDB120', LineStyle='--'
                app.region_Rect.Annotation.LegendInformation.IconDisplayStyle = "off";
                % x=app.Gmax_Theta;
                % y=app.Gmax;
                app.maxGain_XLine = line(app.Axes2, [theta_maxGain theta_maxGain], [app.Rmin maxGain], LineStyle='--', Color='k');
                app.maxGain_YLine = line(app.Axes2, [app.Rmin theta_maxGain], [maxGain maxGain], LineStyle='--', Color='k');
                app.maxGain_XLine.Annotation.LegendInformation.IconDisplayStyle = "off";
                app.maxGain_YLine.Annotation.LegendInformation.IconDisplayStyle = "off";
                
                % text = sprintf('%s\n%s','Line 1','Line 2');
                % label = uilabel('Text',text,'Position',[100 100 100 32]);
                % app.Label_HPBW.Text = ['HPBW:  ' num2str(HPBW) '°'];
                app.Label_HPBW.Text = {['HPBW: ' num2str(HPBW) '°']; ['Max Gain @ ' num2str(theta_maxGain) '°']; ['Max Gain: ' num2str(round(maxGain,2)) 'dBi']};
                % disp(app.Label_HPBW.Text);
                % app.region_Polar = thetaregion(app.pax, [0 2*pi; theta_3dB_1*pi/180 theta_3dB_2*pi/180], 'FaceColor', 'yellow', 'FaceAlpha', 0.3, 'EdgeColor', 'none');
                % app.region_Polar(2).Annotation.LegendInformation.IconDisplayStyle = "off";

                % text(pax, theta_3dB_1*pi/180 + 0.05, Gain_max-3, '-3dB', 'FontSize', 12);
                % text(pax, 0.5, 0.5, sprintf('3dB BW: %.2f degrees', BW), 'HorizontalAlignment', 'center');

                % border lines
                % lower = polarplot(pax, [theta_3dB_1 theta_3dB_1]*pi/180, [min(gain) max_gain], 'k--', 'DisplayName', '3 dB BW');
                % upper = polarplot(pax, [theta_3dB_2 theta_3dB_2]*pi/180, [min(gain) max_gain], 'k--', 'DisplayName', '3 dB BW');
                % 
                % polarplot(pax, [theta_3dB_1 theta_3dB_1]*pi/180, [0 max(gain)], 'b', 'LineWidth', 1.5);
                % polarplot(pax, [theta_3dB_2 theta_3dB_2]*pi/180, [0 max(gain)], 'b', 'LineWidth', 1.5);

                % polarplot(pax, [theta_3dB_1, theta_3dB_1], [0, Gain_3dB], 'g--', 'DisplayName', '3dB BW Start');
                % polarplot(pax, [theta_3dB_2, theta_3dB_2], [0, Gain_3dB], 'g--', 'DisplayName', '3dB BW End');

                % Plot1 = polarplot(app.pax, theta*pi/180, gain,'LineWidth',2,'Color','red',DisplayName='E_Total');
            else
                app.Button_HPBW.BackgroundColor = [0.96,0.96,0.96];
                delete([app.region_Polar, app.region_Rect, app.maxGain_RLine, app.maxGain_XLine, app.maxGain_XLine]);
                app.Label_HPBW.Text = '';
                app.Button_HPBW.Value = 0;
            end
        end

        % Value changed function: Switch
        function SwitchValueChanged(app, event)
            app.plotPattern();
        end

        % Value changed function: DropDown
        function DropDownValueChanged(app, event)
            % value = app.DropDown.Value;
            switch app.DropDown.Value %event.Value
                case {'+X','-X'}
                    app.cutE = app.Theta_Phi_0_180;
                    app.cutH = app.Phi_Theta90;
                    % disp('Selection: X-axis');
                    if app.DropDown.Value == "+X"
                        % xlim(app.Axes2,[0 180]);
                        app.Axes2.XLim = [0 180];
                    else
                        xlim(app.Axes2,[180 360]);
                    end
                case {'+Y','-Y'}
                    % disp('Selection: Y-axis');
                    app.cutE = app.Theta_Phi_90_27;
                    app.cutH = app.Phi_Theta90;
                    if app.DropDown.Value == "+Y"
                        xlim(app.Axes2,[0 180]);
                    else
                        xlim(app.Axes2,[180 360]);
                    end
                otherwise % case app.DropDown.Value == 'Z-axis'
                    % disp('Selection: Z-axis');
                    app.cutE = app.Theta_Phi_0_180;
                    app.cutH = app.Theta_Phi_90_27;
                    xlim(app.Axes2,[0 180]);
            end
            app.Axes2.XTick = min(app.Axes2.XLim):15:max(app.Axes2.XLim);
            % app.Axes2.XTickLabel = string(min(app.Axes2.XLim):15:max(app.Axes2.XLim));
            app.SwitchValueChanged();
        end

        % Button pushed function: Button_ContourPlot
        function Button_ContourPlotPushed(app, event)
            phi_unique = unique(app.Table_DataOut.Data.Phi);
            theta_unique = unique(app.Table_DataOut.Data.Theta);
            [phi_mat, theta_mat] = meshgrid(phi_unique, theta_unique);
            Gain_mat = reshape(app.Table_DataOut.Data.E_Total_dB, length(phi_unique), length(theta_unique));

            epsilon = 0;
            dBmin = floor(min(app.Table_DataOut.Data.E_Total_dB)/10 - epsilon)*10;
            dBmax = ceil(max(app.Table_DataOut.Data.E_Total_dB)/10 + epsilon)*10;

            % fig_Contour = figure('Name','Contour Plot');
            % ax = gca;
            %ax.XLim = [0 10]
            % fig_Contour = uifigure('Name','Contour Plot');
            % ax = uiaxes(fig_Contour);
            % ax.Title.String = 'Antenna Gain Pattern - Contour Plot'; %title(ax,'Antenna Gain Pattern - Contour Plot');
            % ax.XLabel.String = 'Phi (degrees)'; %xlabel(ax,'Phi (degrees)');
            % ax.YLabel.String = 'Theta (degrees)'; %ylabel(ax,'Theta (degrees)');
            % ax.YDir = 'reverse';
            % plotContour = contourf(ax,phi_mat, theta_mat, Gain_mat', 'LineStyle', 'none');
            % colorbar (ax);

            % colorbar;
            figure('Name','Contour Plot','NumberTitle','off');
            plotContour = contourf(phi_mat, theta_mat, Gain_mat', 'LineStyle', 'none');
            clim([dBmin dBmax]); %caxis
            colorbar;
            colormap(jet);
            % fig = gcf;
            % fig.Name = 'Results';
            % fig.NumberTitle = 'off';
            ax = gca ;
            % ax.Title.String = 'Antenna Gain Pattern - Contour Plot'; %title(ax,'Antenna Gain Pattern - Contour Plot');
            ax.XLabel.String = 'Phi (degrees)'; %xlabel(ax,'Phi (degrees)');
            ax.YLabel.String = 'Theta (degrees)'; %ylabel(ax,'Theta (degrees)');
            ax.YDir = 'reverse';
            
        end

        % Callback function
        function SliderLimit_maxValueChanged(app, event)
            app.Rmax = app.SliderLimit_max.Value;
            app.Range.Limits = [app.Rmin app.Rmax];
            app.Range.Value = [app.Rmin app.Rmax];
            MajorTicksStep = 10;
            app.Range.MajorTicks = app.Rmin:MajorTicksStep:app.Rmax;
            app.Range.MajorTickLabels = app.Range.MajorTicks + " dB";
            app.RangeValueChanged();
            app.pax.RTick = (app.Rmin:10:app.Rmax);
            yticks(app.Axes2,(app.Rmin:5:app.Rmax));
        end

        % Callback function
        function SliderLimit_minValueChanged(app, event)
            app.Rmin = app.SliderLimit_min.Value;
            app.Range.Limits = [app.Rmin app.Rmax];
            app.Range.Value = [app.Rmin app.Rmax];
            MajorTicksStep = 10;
            app.Range.MajorTicks = app.Rmin:MajorTicksStep:app.Rmax;
            app.Range.MajorTickLabels = app.Range.MajorTicks + " dB";
            app.RangeValueChanged();
            app.pax.RTick = (app.Rmin:10:app.Rmax);
            yticks(app.Axes2,(app.Rmin:5:app.Rmax));
        end

        % Value changed function: RangeMin
        function RangeMinValueChanged(app, event)
            app.Rmin = app.RangeMin.Value;
            app.Range.Limits = [app.Rmin app.Rmax];
            app.Range.Value = [app.Rmin app.Rmax];
            MajorTicksStep = 10;
            app.Range.MajorTicks = app.Rmin:MajorTicksStep:app.Rmax;
            app.Range.MajorTickLabels = app.Range.MajorTicks + " dB";
            app.RangeValueChanged();
            app.pax.RTick = (app.Rmin:10:app.Rmax);
            yticks(app.Axes2,(app.Rmin:5:app.Rmax));
        end

        % Value changed function: RangeMax
        function RangeMaxValueChanged(app, event)
            app.Rmax = app.RangeMax.Value;
            app.Range.Limits = [app.Rmin app.Rmax];
            app.Range.Value = [app.Rmin app.Rmax];
            MajorTicksStep = 10;
            app.Range.MajorTicks = app.Rmin:MajorTicksStep:app.Rmax;
            app.Range.MajorTickLabels = app.Range.MajorTicks + " dB";
            app.RangeValueChanged();
            app.pax.RTick = (app.Rmin:10:app.Rmax);
            yticks(app.Axes2,(app.Rmin:5:app.Rmax));
        end

        % Button pushed function: Button_ExportCut
        function Button_ExportCutPushed(app, event)
            cd(app.FolderPath);
            % [filename, pathname] = uiputfile('*.csv', 'Save Data as CSV');
            [FileName, PathName, filterindex] = uiputfile({'*.csv', 'CSV Files (*.csv)'; '*.txt', 'Text Files (*.txt)'}, 'Export Cut',[app.InputFileName, '_cut_', app.Switch.Value, '_', char(extractBetween(app.DropDown_step.Value, 'STEP: ', '°')) , '_deg']);
            if isequal(FileName,0) || isequal(PathName,0)
                return; % User canceled the operation
            else
                filepath = fullfile(PathName, FileName);

                if filterindex == 1
                    assignin('base', 'cut2', app.cut);
                    % Export as CSV
                    writetable(app.cut, filepath);
                elseif filterindex == 2
                    % Export as TXT
                    writetable(app.cut, filepath, 'Delimiter', '\t');
                end
            end
        end

        % Value changed function: CheckBox_Et
        function CheckBox_EtValueChanged(app, event)
            % Control E_Total
            if app.CheckBox_Et.Value
                app.Plot1.Visible = 'on';
            else
                app.Plot1.Visible = 'off';
            end
        end

        % Value changed function: CheckBox_Er
        function CheckBox_ErValueChanged(app, event)
            % Control E_RCP
            if app.CheckBox_Er.Value
                app.Plot1_2.Visible = 'on';
            else
                app.Plot1_2.Visible = 'off';
            end
        end

        % Value changed function: CheckBox_El
        function CheckBox_ElValueChanged(app, event)
            % Control E_LCP
            if app.CheckBox_El.Value
                app.Plot1_3.Visible = 'on';
            else
                app.Plot1_3.Visible = 'off';
            end
        end

        % Value changed function: DropDown_step
        function DropDown_stepValueChanged(app, event)
            app.processData();
            assignin('base', 'cut', app.cut);
        end

        % Value changed function: Format
        function FormatValueChanged(app, event)
            if ~isempty(app.FilePath) % or if ~isempty(app.Table_DataIn.Data)
                app.processData();
            end
        end
    end

    % Component initialization
    methods (Access = private)

        % Create UIFigure and components
        function createComponents(app)

            % Create UIFigure and hide until all components are created
            app.UIFigure = uifigure('Visible', 'off');
            app.UIFigure.Position = [100 100 1150 750];
            app.UIFigure.Name = 'PatternReader';
            app.UIFigure.WindowState = 'maximized';

            % Create Grid
            app.Grid = uigridlayout(app.UIFigure);
            app.Grid.ColumnWidth = {'1x', '1x', '1x', '1x', '1x', '1x', '1x', '1x'};
            app.Grid.RowHeight = {'fit', 'fit', 'fit', '1x', '1x', 'fit', '1x', '1x', 'fit'};

            % Create Button_Load
            app.Button_Load = uibutton(app.Grid, 'push');
            app.Button_Load.ButtonPushedFcn = createCallbackFcn(app, @Button_LoadPushed, true);
            app.Button_Load.FontSize = 14;
            app.Button_Load.FontWeight = 'bold';
            app.Button_Load.Layout.Row = 2;
            app.Button_Load.Layout.Column = [3 7];
            app.Button_Load.Text = 'Load Antenna Pattern';

            % Create EditField_Loss
            app.EditField_Loss = uispinner(app.Grid);
            app.EditField_Loss.Step = 0.1;
            app.EditField_Loss.Limits = [0 Inf];
            app.EditField_Loss.ValueChangedFcn = createCallbackFcn(app, @EditField_LossValueChanged2, true);
            app.EditField_Loss.Visible = 'off';
            app.EditField_Loss.Layout.Row = 3;
            app.EditField_Loss.Layout.Column = 4;
            app.EditField_Loss.Value = 0.3;

            % Create Panel_Polar
            app.Panel_Polar = uipanel(app.Grid);
            app.Panel_Polar.TitlePosition = 'centertop';
            app.Panel_Polar.Title = 'Polar Plot';
            app.Panel_Polar.Visible = 'off';
            app.Panel_Polar.BackgroundColor = [0.9412 0.9412 0.9412];
            app.Panel_Polar.Layout.Row = [4 5];
            app.Panel_Polar.Layout.Column = [1 4];
            app.Panel_Polar.FontWeight = 'bold';

            % Create Grid_Polar
            app.Grid_Polar = uigridlayout(app.Panel_Polar);
            app.Grid_Polar.ColumnWidth = {'fit', '0.26x', '1x', '0.23x'};
            app.Grid_Polar.RowHeight = {'fit', '0.25x', '1x', 'fit'};

            % Create Range
            app.Range = uislider(app.Grid_Polar, 'range');
            app.Range.Orientation = 'vertical';
            app.Range.ValueChangedFcn = createCallbackFcn(app, @RangeValueChanged, true);
            app.Range.Step = 1;
            app.Range.Layout.Row = [2 3];
            app.Range.Layout.Column = 1;

            % Create Button_HPBW
            app.Button_HPBW = uibutton(app.Grid_Polar, 'state');
            app.Button_HPBW.ValueChangedFcn = createCallbackFcn(app, @Button_HPBWValueChanged, true);
            app.Button_HPBW.IconAlignment = 'center';
            app.Button_HPBW.Text = 'HPBW';
            app.Button_HPBW.FontWeight = 'bold';
            app.Button_HPBW.Layout.Row = 1;
            app.Button_HPBW.Layout.Column = 4;

            % Create Label_HPBW
            app.Label_HPBW = uilabel(app.Grid_Polar);
            app.Label_HPBW.HorizontalAlignment = 'center';
            app.Label_HPBW.FontWeight = 'bold';
            app.Label_HPBW.Layout.Row = 2;
            app.Label_HPBW.Layout.Column = 4;
            app.Label_HPBW.Text = '';

            % Create RangeMin
            app.RangeMin = uispinner(app.Grid_Polar);
            app.RangeMin.Step = 5;
            app.RangeMin.ValueChangedFcn = createCallbackFcn(app, @RangeMinValueChanged, true);
            app.RangeMin.Layout.Row = 4;
            app.RangeMin.Layout.Column = 1;

            % Create RangeMax
            app.RangeMax = uispinner(app.Grid_Polar);
            app.RangeMax.Step = 5;
            app.RangeMax.ValueChangedFcn = createCallbackFcn(app, @RangeMaxValueChanged, true);
            app.RangeMax.Layout.Row = 1;
            app.RangeMax.Layout.Column = 1;

            % Create Button_ExportCut
            app.Button_ExportCut = uibutton(app.Grid_Polar, 'push');
            app.Button_ExportCut.ButtonPushedFcn = createCallbackFcn(app, @Button_ExportCutPushed, true);
            app.Button_ExportCut.FontWeight = 'bold';
            app.Button_ExportCut.Layout.Row = 4;
            app.Button_ExportCut.Layout.Column = 4;
            app.Button_ExportCut.Text = 'Export Cut';

            % Create GridLayout
            app.GridLayout = uigridlayout(app.Grid_Polar);
            app.GridLayout.ColumnWidth = {'1x'};
            app.GridLayout.RowHeight = {'1x', '1x', '1x'};
            app.GridLayout.Layout.Row = 3;
            app.GridLayout.Layout.Column = 4;

            % Create CheckBox_El
            app.CheckBox_El = uicheckbox(app.GridLayout);
            app.CheckBox_El.ValueChangedFcn = createCallbackFcn(app, @CheckBox_ElValueChanged, true);
            app.CheckBox_El.Text = 'E_LCP';
            app.CheckBox_El.Layout.Row = 3;
            app.CheckBox_El.Layout.Column = 1;
            app.CheckBox_El.Value = true;

            % Create CheckBox_Er
            app.CheckBox_Er = uicheckbox(app.GridLayout);
            app.CheckBox_Er.ValueChangedFcn = createCallbackFcn(app, @CheckBox_ErValueChanged, true);
            app.CheckBox_Er.Text = 'E_RCP';
            app.CheckBox_Er.Layout.Row = 2;
            app.CheckBox_Er.Layout.Column = 1;
            app.CheckBox_Er.Value = true;

            % Create CheckBox_Et
            app.CheckBox_Et = uicheckbox(app.GridLayout);
            app.CheckBox_Et.ValueChangedFcn = createCallbackFcn(app, @CheckBox_EtValueChanged, true);
            app.CheckBox_Et.Text = 'E_Total';
            app.CheckBox_Et.Layout.Row = 1;
            app.CheckBox_Et.Layout.Column = 1;
            app.CheckBox_Et.Value = true;

            % Create EditField_Rw
            app.EditField_Rw = uispinner(app.Grid);
            app.EditField_Rw.ValueChangedFcn = createCallbackFcn(app, @EditField_RwValueChanged2, true);
            app.EditField_Rw.Visible = 'off';
            app.EditField_Rw.Layout.Row = 3;
            app.EditField_Rw.Layout.Column = 7;
            app.EditField_Rw.Value = 6;

            % Create Panel_Rect
            app.Panel_Rect = uipanel(app.Grid);
            app.Panel_Rect.TitlePosition = 'centertop';
            app.Panel_Rect.Title = 'Rectangular Plot';
            app.Panel_Rect.Visible = 'off';
            app.Panel_Rect.BackgroundColor = [0.9412 0.9412 0.9412];
            app.Panel_Rect.Layout.Row = [4 5];
            app.Panel_Rect.Layout.Column = [5 8];
            app.Panel_Rect.FontWeight = 'bold';

            % Create Grid_Rect
            app.Grid_Rect = uigridlayout(app.Panel_Rect);

            % Create Axes2
            app.Axes2 = uiaxes(app.Grid_Rect);
            xlabel(app.Axes2, 'Theta (degree)')
            ylabel(app.Axes2, 'Magnitude (dB)')
            zlabel(app.Axes2, 'Z')
            app.Axes2.XLim = [0 180];
            app.Axes2.XTick = [0 15 30 45 60 75 90 105 120 135 150 165 180];
            app.Axes2.Box = 'on';
            app.Axes2.Layout.Row = [1 2];
            app.Axes2.Layout.Column = [1 2];
            app.Axes2.Visible = 'off';

            % Create Label_max
            app.Label_max = uilabel(app.Grid);
            app.Label_max.BackgroundColor = [0.8 0.8 0.8];
            app.Label_max.HorizontalAlignment = 'center';
            app.Label_max.FontSize = 13;
            app.Label_max.FontWeight = 'bold';
            app.Label_max.Visible = 'off';
            app.Label_max.Layout.Row = 6;
            app.Label_max.Layout.Column = [1 2];

            % Create Button_3dPlot
            app.Button_3dPlot = uibutton(app.Grid, 'push');
            app.Button_3dPlot.ButtonPushedFcn = createCallbackFcn(app, @Button_3dPlotPushed, true);
            app.Button_3dPlot.FontWeight = 'bold';
            app.Button_3dPlot.Visible = 'off';
            app.Button_3dPlot.Layout.Row = 6;
            app.Button_3dPlot.Layout.Column = [4 5];
            app.Button_3dPlot.Text = 'Plot3D';

            % Create Button_ContourPlot
            app.Button_ContourPlot = uibutton(app.Grid, 'push');
            app.Button_ContourPlot.ButtonPushedFcn = createCallbackFcn(app, @Button_ContourPlotPushed, true);
            app.Button_ContourPlot.FontWeight = 'bold';
            app.Button_ContourPlot.Visible = 'off';
            app.Button_ContourPlot.Layout.Row = 6;
            app.Button_ContourPlot.Layout.Column = 3;
            app.Button_ContourPlot.Text = 'Contour Plot';

            % Create Label_Pol
            app.Label_Pol = uilabel(app.Grid);
            app.Label_Pol.BackgroundColor = [0.8 0.8 0.8];
            app.Label_Pol.HorizontalAlignment = 'center';
            app.Label_Pol.FontSize = 13;
            app.Label_Pol.FontWeight = 'bold';
            app.Label_Pol.Visible = 'off';
            app.Label_Pol.Layout.Row = 6;
            app.Label_Pol.Layout.Column = [7 8];
            app.Label_Pol.Text = '';

            % Create Table_DataIn
            app.Table_DataIn = uitable(app.Grid);
            app.Table_DataIn.BackgroundColor = [1 1 1];
            app.Table_DataIn.ColumnName = {'Theta'; 'Phi'; 'E-TH-DB'; 'E-PH-DB'; 'E-TH-DG'; 'E-PH-DG'};
            app.Table_DataIn.ColumnRearrangeable = 'on';
            app.Table_DataIn.RowName = {};
            app.Table_DataIn.ColumnSortable = true;
            app.Table_DataIn.Enable = 'off';
            app.Table_DataIn.Visible = 'off';
            app.Table_DataIn.Layout.Row = [7 8];
            app.Table_DataIn.Layout.Column = [1 4];

            % Create Table_DataOut
            app.Table_DataOut = uitable(app.Grid);
            app.Table_DataOut.BackgroundColor = [1 1 1];
            app.Table_DataOut.ColumnName = {'Theta'; 'Phi'; 'E_Total_dB'; 'E_RCP_dB'; 'E_LCP_dB'; 'AR_dB'; 'PLF_dB'; 'E_Polarized'};
            app.Table_DataOut.ColumnRearrangeable = 'on';
            app.Table_DataOut.RowName = {};
            app.Table_DataOut.ColumnSortable = true;
            app.Table_DataOut.Visible = 'off';
            app.Table_DataOut.Layout.Row = [7 8];
            app.Table_DataOut.Layout.Column = [5 8];

            % Create Export_Input
            app.Export_Input = uibutton(app.Grid, 'push');
            app.Export_Input.ButtonPushedFcn = createCallbackFcn(app, @Export_InputButtonPushed, true);
            app.Export_Input.FontWeight = 'bold';
            app.Export_Input.Visible = 'off';
            app.Export_Input.Layout.Row = 9;
            app.Export_Input.Layout.Column = [2 3];
            app.Export_Input.Text = 'Export Theta-Phi Antenna Pattern';

            % Create Export_Output
            app.Export_Output = uibutton(app.Grid, 'push');
            app.Export_Output.ButtonPushedFcn = createCallbackFcn(app, @Export_OutputButtonPushed, true);
            app.Export_Output.FontWeight = 'bold';
            app.Export_Output.Visible = 'off';
            app.Export_Output.Layout.Row = 9;
            app.Export_Output.Layout.Column = [6 7];
            app.Export_Output.Text = 'Export Gain Antenna Pattern';

            % Create Switch
            app.Switch = uiswitch(app.Grid, 'slider');
            app.Switch.Items = {'E-Plane', 'H-Plane'};
            app.Switch.ValueChangedFcn = createCallbackFcn(app, @SwitchValueChanged, true);
            app.Switch.Visible = 'off';
            app.Switch.Layout.Row = 3;
            app.Switch.Layout.Column = 1;
            app.Switch.Value = 'E-Plane';

            % Create DropDown
            app.DropDown = uidropdown(app.Grid);
            app.DropDown.Items = {'Orientation:', '+Z', '-Z', '+X', '-X', '+Y', '-Y'};
            app.DropDown.ValueChangedFcn = createCallbackFcn(app, @DropDownValueChanged, true);
            app.DropDown.Visible = 'off';
            app.DropDown.Layout.Row = 3;
            app.DropDown.Layout.Column = 2;
            app.DropDown.Value = 'Orientation:';

            % Create Button_Coverage
            app.Button_Coverage = uibutton(app.Grid, 'push');
            app.Button_Coverage.FontWeight = 'bold';
            app.Button_Coverage.Visible = 'off';
            app.Button_Coverage.Layout.Row = 6;
            app.Button_Coverage.Layout.Column = 6;
            app.Button_Coverage.Text = 'Coverage';

            % Create LossindBLabel
            app.LossindBLabel = uilabel(app.Grid);
            app.LossindBLabel.HorizontalAlignment = 'right';
            app.LossindBLabel.Visible = 'off';
            app.LossindBLabel.Layout.Row = 3;
            app.LossindBLabel.Layout.Column = 3;
            app.LossindBLabel.Text = 'Loss (in dB)';

            % Create IncidentWaveARRwPLFLabel
            app.IncidentWaveARRwPLFLabel = uilabel(app.Grid);
            app.IncidentWaveARRwPLFLabel.HorizontalAlignment = 'right';
            app.IncidentWaveARRwPLFLabel.Visible = 'off';
            app.IncidentWaveARRwPLFLabel.Layout.Row = 3;
            app.IncidentWaveARRwPLFLabel.Layout.Column = [5 6];
            app.IncidentWaveARRwPLFLabel.Text = '''Incident Wave AxialRatio (''Rw in dB) [for PLF]''';

            % Create FormatLabel
            app.FormatLabel = uilabel(app.Grid);
            app.FormatLabel.HorizontalAlignment = 'right';
            app.FormatLabel.Layout.Row = 2;
            app.FormatLabel.Layout.Column = 1;
            app.FormatLabel.Text = 'Format:';

            % Create Format
            app.Format = uidropdown(app.Grid);
            app.Format.Items = {'Eθ, Eϕ [mag_dB,phase°]', 'Eθ, Eϕ [Re,Im]', 'Erh, Elh [mag,phase°]', 'Erh, Elh [Re,Im]'};
            app.Format.ItemsData = [1 2 3 4];
            app.Format.ValueChangedFcn = createCallbackFcn(app, @FormatValueChanged, true);
            app.Format.Layout.Row = 2;
            app.Format.Layout.Column = 2;
            app.Format.Value = 1;

            % Create DropDown_step
            app.DropDown_step = uidropdown(app.Grid);
            app.DropDown_step.Items = {'STEP', 'STEP: 1'};
            app.DropDown_step.ValueChangedFcn = createCallbackFcn(app, @DropDown_stepValueChanged, true);
            app.DropDown_step.Enable = 'off';
            app.DropDown_step.Visible = 'off';
            app.DropDown_step.Placeholder = 'STEP';
            app.DropDown_step.Layout.Row = 9;
            app.DropDown_step.Layout.Column = [4 5];
            app.DropDown_step.Value = 'STEP';

            % Show the figure after all components are created
            app.UIFigure.Visible = 'on';
        end
    end

    % App creation and deletion
    methods (Access = public)

        % Construct app
        function app = Antenna_Pattern_Ppocessing_2

            % Create UIFigure and components
            createComponents(app)

            % Register the app with App Designer
            registerApp(app, app.UIFigure)

            % Execute the startup function
            runStartupFcn(app, @startupFcn)

            if nargout == 0
                clear app
            end
        end

        % Code that executes before app deletion
        function delete(app)

            % Delete UIFigure when app is deleted
            delete(app.UIFigure)
        end
    end
end