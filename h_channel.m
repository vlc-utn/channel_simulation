function [h, delay] = h_channel(r_s, n_s, m, h_s, r_r, n_r, A, FOV)
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
    %
    % Outputs:
    %   - h = Attenuation of the optical channel.
    %   - delay = Temporal delay to travel from the sender to the receiver.
    arguments(Input)
        r_s (:, 3) double
        n_s (:, 3) double
        m double
        h_s (:,1) double
        r_r (:, 3) double
        n_r (:, 3) double
        A double
        FOV double
    end
    arguments(Output)
        h (:,1) double
        delay (:,1) double
    end

    h = zeros(height(r_r), 1);

    for j=1:1:height(r_s)
        % Pre-allocate vectors
        distance = zeros(height(r_r), 1);
        cos_emitter = zeros(size(distance));
        cos_receiver = zeros(size(distance));
    
        % Vector operations
        for i=1:1:length(r_r)
            distance(i) = norm(r_s(j,:) - r_r(i,:));
            cos_emitter(i) = dot(n_s(j,:), (r_r(i,:) - r_s(j,:)) ./ distance(i));
            cos_receiver(i) = dot(n_r(i,:), (r_s(j,:) - r_r(i,:)) ./ distance(i));
        end
        cos_emitter(cos_emitter < 0) = 0;

        h = h + h_s(j,:) .* ((m+1) / (2*pi)) .* cos_emitter.^m .* A .* cos_receiver ...
            .* rect(acosd(cos_receiver) / FOV) ./ (distance.^2);
    end
    
    

    delay = distance ./ physconst("LightSpeed");

end