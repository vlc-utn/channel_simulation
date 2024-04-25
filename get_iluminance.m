function [Iluminance] = get_iluminance(r_s, n_s, m, r_r, n_r, I0)
    %H_CHANNEL. Impulse response from the optical channel.
    %
    % Args:
    %   - r_s = (x,y,z) position of the senders.
    %   - n_s = (x,y,z) orientation of the senders.
    %   - m = Lambert mode of the sender.
    %   - r_r = (x,y,z) position of the receivers.
    %   - n_r = (x,y,z) position of the receivers.
    %   - I0 = [lm] Total luminic power, or normalized luminic intensity. 
    %
    % Outputs:
    %   - Iluminance = [lx] Iluminance.
    arguments(Input)
        r_s (:, 3) double
        n_s (:, 3) double
        m double
        r_r (:, 3) double
        n_r (:, 3) double
        I0 double
    end
    arguments(Output)
        Iluminance (:,:) double
    end

    Iluminance = zeros(height(r_r), 1);

    for j=1:1:height(r_s)
        % Vector operations
        distance = vecnorm(r_s(j,:) - r_r, 2, 2);
        distance(distance == 0) = 1;
        cos_emitter = dot(repmat(n_s(j,:), height(r_r), 1), (r_r - r_s(j,:)) ./ distance, 2);
        cos_receiver = dot(n_r, (r_s(j,:) - r_r) ./ distance, 2);
        cos_emitter(cos_emitter < 0) = 0;

        % TODO agregar coseno del transmisor
        Iluminance = Iluminance + ((m+1) / (2*pi)) .* I0 .* cos_emitter.^m ...
            .* cos_receiver ./ (distance.^2);
    end
end