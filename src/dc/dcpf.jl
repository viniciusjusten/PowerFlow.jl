function run_dcpf(input)
    n = num_buses(input)
    B = dc_susceptance_matrix(input)

    slack_idx = findfirst(i -> bus_type(input, i) == BusType.Slack, 1:n)
    non_slack = [i for i in 1:n if i != slack_idx]

    P = [bus_active_power_injection(input, i) for i in non_slack]
    B_reduced = B[non_slack, non_slack]
    angle_reduced = B_reduced \ P

    angle = zeros(n)
    angle[slack_idx] = bus_voltage_angle(input, slack_idx)
    for (k, i) in enumerate(non_slack)
        angle[i] = angle_reduced[k]
    end

    Pg = [generator_active_power(input, g) for g in 1:num_generators(input)]

    return DCSolution(angle, Pg, NaN, nothing)
end
