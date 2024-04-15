function [h, delay] = h_channel(r_s, n_s, m, r_r, n_r, A, FOV)
    %H_CHANNEL. Impulse response from the optical channel.
    %
    % Args:
    %   - r_s = (x,y,z) position of the sender.
    %   - n_s = (x,y,z) orientation of the sender.
    %   - m = Lambert mode of the sender.
    %   - r_r = (x,y,z) position of the receivers.
    %   - n_r = (x,y,z) position of the receivers.
    %   - A = area of all the receivers.
    %   - FOV = Field of view of all the receivers. Incident beams with an
    %   angle greater than the FOV will be discarded.
    %
    % Outputs:
    %   - h = Attenuation of the optical channel.
    %   - delay = Temporal delay to travel from the sender to the receiver.
    arguments(Input)
        r_s (1, 3) double
        n_s (1, 3) double
        m double
        r_r (:, 3) double
        n_r (:, 3) double
        A double
        FOV double
    end
    arguments(Output)
        h (:,1) double
        delay (:,1) double
    end

    % Pre-allocate vectors
    distance = zeros(length(r_r), 1);
    cos_emitter = zeros(size(distance));
    cos_receiver = zeros(size(distance));

    % Vector operations
    for i=1:1:length(r_r)
        distance(i) = norm(r_s - r_r(i,:));
        cos_emitter(i) = dot(n_s, (r_r(i,:) - r_s) ./ distance(i));
        cos_receiver(i) = dot(n_r(i,:), (r_s - r_r(i,:)) ./ distance(i));
    end

    % LOS channel response
    h = ((m+1) / (2*pi)) .* cos_emitter.^m .* A .* cos_receiver ...
        .* rect(acosd(cos_receiver) / FOV) ./ (distance.^2);

    delay = distance ./ physconst("LightSpeed");

end