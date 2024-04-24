function [h] = h_channel(r_s, n_s, m, h_s, r_r, n_r, A, FOV, t_vector)
    %H_CHANNEL. Impulse response from the optical channel.
    %
    % Args:
    %   - r_s = (x,y,z) position of the senders.
    %   - n_s = (x,y,z) orientation of the senders.
    %   - m = Lambert mode of the sender.
    %   - h_s = Impulse response carried from the previous sender.
    %   - r_r = (x,y,z) position of the receivers.
    %   - n_r = (x,y,z) position of the receivers.
    %   - A = area of all the receivers.
    %   - FOV = Field of view of all the receivers. Incident beams with an
    %   angle greater than the FOV will be discarded.
    %   - t_vector: Temporal vector.
    %
    % Outputs:
    %   - h = Attenuation of the optical channel.
    %   - delay = Temporal delay to travel from the sender to the receiver.
    arguments(Input)
        r_s (:, 3) double
        n_s (:, 3) double
        m double
        h_s (:, :) double
        r_r (:, 3) double
        n_r (:, 3) double
        A double
        FOV double
        t_vector (1, :) double
    end
    arguments(Output)
        h (:,:) double
    end

    h = zeros(height(r_r), length(t_vector));

    for j=1:1:height(r_s)
        % Pre-allocate vectors
        h_aux = zeros(height(r_r), length(t_vector));

        % Vector operations
        distance = vecnorm(r_s(j,:) - r_r, 2, 2);
        distance(distance == 0) = 1;
        cos_emitter = dot(repmat(n_s(j,:), height(r_r), 1), (r_r - r_s(j,:)) ./ distance, 2);
        cos_receiver = dot(n_r, (r_s(j,:) - r_r) ./ distance, 2);
        cos_emitter(cos_emitter < 0) = 0;

        channel_att = ((m+1) / (2*pi)) .* cos_emitter.^m .* A .* cos_receiver ...
            .* rect(acosd(cos_receiver) / FOV) ./ (distance.^2);

        [~, index] = min(abs(distance/physconst("LightSpeed") - t_vector), [], 2);

        for i=1:1:length(r_r)
            h_aux(i,index(i)) = channel_att(i);
        end

        h_conv = conv2(h_aux, h_s(j,:), "full");
        h_aux = h_conv(:, 1:length(t_vector));
        h = h + h_aux;
    end
end