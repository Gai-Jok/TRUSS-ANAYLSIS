function truss_analysis_enhanced()
    clc; clear; close all;

    %% Input Parameters
    disp('=== TRUSS ANALYSIS PROGRAM ===');

    % Node information
    num_nodes = input('Enter the number of nodes: ');
    nodes = zeros(num_nodes, 2);
    disp('Enter node coordinates as [x y]:');
    for i = 1:num_nodes
        nodes(i,:) = input(sprintf('  Node %d: ', i));
    end

    % Element information
    num_elements = input('Enter the number of elements: ');
    elements = zeros(num_elements, 2);
    disp('Enter element connectivity as [Node1 Node2]:');
    for i = 1:num_elements
        elements(i,:) = input(sprintf('  Element %d: ', i));
    end

    % Material properties
    num_E = input('How many different Young''s moduli are there? ');
    E_values = zeros(num_E, 1);
    for i = 1:num_E
        E_values(i) = input(sprintf('  Young''s Modulus %d (Pa): ', i));
    end

    num_A = input('How many different cross-sectional areas are there? ');
    A_values = zeros(num_A, 1);
    for i = 1:num_A
        A_values(i) = input(sprintf('  Cross-sectional area %d (m^2): ', i));
    end

    % Assign properties to elements
    E = zeros(num_elements, 1);
    A = zeros(num_elements, 1);
    disp('Assign properties to each element:');
    for i = 1:num_elements
        fprintf('\nElement %d (Nodes %d-%d):\n', i, elements(i,1), elements(i,2));

        % Young's modulus assignment
        if num_E > 1
            fprintf('Available Young''s moduli:\n');
            for j = 1:num_E
                fprintf('  %d: %.2e Pa\n', j, E_values(j));
            end
            e_idx = input('  Select Young''s modulus (enter index): ');
            E(i) = E_values(e_idx);
        else
            E(i) = E_values(1);
        end

        % Area assignment
        if num_A > 1
            fprintf('Available cross-sectional areas:\n');
            for j = 1:num_A
                fprintf('  %d: %.2e m²\n', j, A_values(j));
            end
            a_idx = input('  Select cross-sectional area (enter index): ');
            A(i) = A_values(a_idx);
        else
            A(i) = A_values(1);
        end
    end

    % Forces
    forces = zeros(2*num_nodes, 1);
    disp('\nEnter forces at each node as [Fx Fy] (enter 0 if no force):');
    for i = 1:num_nodes
        f = input(sprintf('  Force at Node %d: ', i));
        forces(2*i-1:2*i) = f(:);
    end

    % Boundary conditions
    fixed_dofs = zeros(2*num_nodes, 1);
    disp('\nEnter boundary conditions as [Ux Uy] (1 = fixed, 0 = free):');
    for i = 1:num_nodes
        bc = input(sprintf('  Node %d (Fixed/Free): ', i));
        fixed_dofs(2*i-1:2*i) = bc(:);
    end

    %% Finite Element Analysis
    % Global stiffness matrix
    K = zeros(2*num_nodes, 2*num_nodes);

    for e = 1:num_elements
        n1 = elements(e,1);
        n2 = elements(e,2);
        x1 = nodes(n1,1); y1 = nodes(n1,2);
        x2 = nodes(n2,1); y2 = nodes(n2,2);

        L = sqrt((x2-x1)^2 + (y2-y1)^2);
        c = (x2-x1)/L; s = (y2-y1)/L;

        % Element stiffness matrix
        ke = (E(e)*A(e)/L) * [c^2 c*s -c^2 -c*s;
                             c*s s^2 -c*s -s^2;
                            -c^2 -c*s c^2 c*s;
                            -c*s -s^2 c*s s^2];

        % Assembly
        dof_map = [2*n1-1, 2*n1, 2*n2-1, 2*n2];
        K(dof_map, dof_map) = K(dof_map, dof_map) + ke;
    end

    % Apply boundary conditions
    free_dofs = find(~fixed_dofs);
    fixed_dofs_indices = find(fixed_dofs); % Renamed to avoid confusion
    K_red = K(free_dofs, free_dofs);
    F_red = forces(free_dofs);

    % Solve for displacements
    U_red = K_red \ F_red;

    % Full displacement vector
    U = zeros(2*num_nodes, 1);
    U(free_dofs) = U_red;

    % Calculate reactions
    R = K * U - forces; % Correct reaction calculation (R = KU - F_external)
    reactions = R(fixed_dofs_indices);

    % Member forces and stresses
    forces_member = zeros(num_elements, 1);
    stresses = zeros(num_elements, 1);

    for e = 1:num_elements
        n1 = elements(e,1);
        n2 = elements(e,2);
        x1 = nodes(n1,1); y1 = nodes(n1,2);
        x2 = nodes(n2,1); y2 = nodes(n2,2);

        L = sqrt((x2-x1)^2 + (y2-y1)^2);
        c = (x2-x1)/L; s = (y2-y1)/L;

        ue = [U(2*n1-1); U(2*n1); U(2*n2-1); U(2*n2)];
        forces_member(e) = (E(e)*A(e)/L) * [-c -s c s] * ue;
        stresses(e) = forces_member(e)/A(e);
    end

    %% Results Display
    disp(' ');
    disp('=== ANALYSIS RESULTS ===');

    % Nodal displacements
    disp('Nodal Displacements (m):');
    for i = 1:num_nodes
        fprintf('Node %d: Ux = %.4e, Uy = %.4e\n', i, U(2*i-1), U(2*i));
    end

    % Reactions
    disp(' ');
    disp('Reactions at fixed nodes (N):');
    for i = 1:length(fixed_dofs_indices)
        dof_idx = fixed_dofs_indices(i);
        node = ceil(dof_idx/2);
        dof_type = mod(dof_idx,2);
        if dof_type == 1
            fprintf('Node %d, Rx = %.2f N\n', node, reactions(i));
        else
            fprintf('Node %d, Ry = %.2f N\n', node, reactions(i));
        end
    end

    % Member forces
    disp(' ');
    disp('Member Forces (N) and Stresses (Pa):');
    for i = 1:num_elements
        fprintf('Element %d (Nodes %d-%d): Force = %.2f N, Stress = %.2e Pa\n',...
                i, elements(i,1), elements(i,2), forces_member(i), stresses(i));
    end

    %% Visualization (Enhanced for Labels and Arrow Scaling)
    figure('Name', 'Truss Analysis Results', 'Color', 'w');
    hold on; axis equal; grid on; box on;

    % --- SCALING CALCULATIONS ---
    % Determine truss size to scale graphical elements relative to the structure
    x_range = range(nodes(:,1));
    y_range = range(nodes(:,2));
    truss_size = max([x_range, y_range]);
    if truss_size == 0, truss_size = 1; end

    % Displacement scaling (exaggeration)
    max_disp = max(abs(U));
    if max_disp == 0
        disp_scale = 1;
    else
        disp_scale = 0.15 * truss_size / max_disp; % Deformed shape moves ~15% of truss size
    end

    % Arrow scaling: Make largest reaction arrow ~15% of truss size
    max_react = max(abs(reactions));
    if max_react < 1e-5, max_react = 1; end % Prevent divide by zero
    arrow_len_factor = (0.15 * truss_size) / max_react;

    % --- 1. PLOT UNDEFORMED SHAPE ---
    for i = 1:num_elements
        n1 = elements(i,1); n2 = elements(i,2);
        h_orig = plot([nodes(n1,1), nodes(n2,1)], [nodes(n1,2), nodes(n2,2)], ...
             'k:', 'LineWidth', 1, 'DisplayName', 'Undeformed');
    end

    % --- 2. PLOT DEFORMED SHAPE ---
    deformed_nodes = nodes + disp_scale*reshape(U,2,num_nodes)';
    for i = 1:num_elements
        n1 = elements(i,1); n2 = elements(i,2);
        h_def = plot([deformed_nodes(n1,1), deformed_nodes(n2,1)], ...
                     [deformed_nodes(n1,2), deformed_nodes(n2,2)], ...
                     'r-', 'LineWidth', 2.5, 'DisplayName', 'Deformed');
    end

    % --- 3. PLOT NODES ---
    plot(nodes(:,1), nodes(:,2), 'ko', 'MarkerFaceColor', 'w', 'MarkerSize', 6, 'HandleVisibility','off');

    % --- 4. TENSION / COMPRESSION MARKERS ---
    h_tens = []; h_comp = [];
    for i = 1:num_elements
        mid = mean(nodes(elements(i,:),:));
        if forces_member(i) > 1e-3
            h_tens = plot(mid(1), mid(2), 'o', 'MarkerSize', 10, ...
                'MarkerFaceColor', [1 0.6 0.6], 'MarkerEdgeColor', 'r', 'DisplayName', 'Tension');
        elseif forces_member(i) < -1e-3
            h_comp = plot(mid(1), mid(2), 'o', 'MarkerSize', 10, ...
                'MarkerFaceColor', [0.6 0.6 1], 'MarkerEdgeColor', 'b', 'DisplayName', 'Compression');
        end
    end

    % --- 5. REACTION ARROWS & LABELS ---
    h_react = [];
    for i = 1:length(fixed_dofs_indices)
        val = reactions(i);
        if abs(val) < 0.1, continue; end % Skip negligible reactions

        dof = fixed_dofs_indices(i);
        node_idx = ceil(dof/2);
        x = nodes(node_idx, 1);
        y = nodes(node_idx, 2);

        % Direction Vector
        if mod(dof, 2) == 1 % X-Force
            u_dir = val; v_dir = 0;
        else % Y-Force
            u_dir = 0; v_dir = val;
        end

        % Normalize magnitude for visualization, keep sign
        draw_length = abs(val) * arrow_len_factor;
        u_draw = sign(u_dir) * draw_length;
        v_draw = sign(v_dir) * draw_length;

        % Plot Arrow (Blue, Thick)
        % 'AutoScale','off' ensures we control exact length
        h_react = quiver(x, y, u_draw, v_draw, 0, ...
            'Color', [0 0 0.8], 'LineWidth', 2.5, 'MaxHeadSize', 0.5, ...
            'DisplayName', 'Reactions');

        % Plot Text Label
        label_txt = sprintf('%.1f N', val);

        % Offset text slightly past the arrow tip
        txt_x = x + u_draw * 1.1;
        txt_y = y + v_draw * 1.1;

        % Alignment logic based on direction
        halign = 'center'; valign = 'middle';
        if abs(u_draw) > abs(v_draw) % Horizontal arrow
            if u_draw > 0, halign = 'left'; else, halign = 'right'; end
        else % Vertical arrow
            if v_draw > 0, valign = 'bottom'; else, valign = 'top'; end
        end

        text(txt_x, txt_y, label_txt, ...
            'Color', 'b', 'FontWeight', 'bold', 'FontSize', 9, ...
            'HorizontalAlignment', halign, 'VerticalAlignment', valign, ...
            'BackgroundColor', 'w', 'EdgeColor', 'none', 'Margin', 1);
    end

    % --- 6. LEGEND & FORMATTING ---
    % Collect handles for legend to avoid duplicates
    legend_handles = [h_orig(1), h_def(1)];
    if ~isempty(h_tens), legend_handles(end+1) = h_tens(1); end
    if ~isempty(h_comp), legend_handles(end+1) = h_comp(1); end
    if ~isempty(h_react), legend_handles(end+1) = h_react(1); end

    legend(legend_handles, 'Location', 'bestoutside');

    % Pad axis limits so arrows/text don't get cut off
    padding = 0.2 * truss_size;
    xlim([min(nodes(:,1))-padding, max(nodes(:,1))+padding]);
    ylim([min(nodes(:,2))-padding, max(nodes(:,2))+padding]);

    title('Truss Analysis: Deformed Shape & Reactions');
    xlabel('X Coordinate (m)');
    ylabel('Y Coordinate (m)');

    hold off;
end
