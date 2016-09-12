classdef C < handle
    properties (Constant)
        % Gesture movement
        G_L = 1;
        G_R = 2;
        G_U = 3;
        G_D = 4;
        G_NONE = 5;
        
        % Wearing wrists
        W_L = 1;
        W_R = 2;
    end

    methods (Static)
        function type_str = get_gesture_movement_name(label)
            if label == C.G_L
                type_str = 'Moving left';
            elseif label == C.G_R
                type_str = 'Moving right';
            elseif label == C.G_U
                type_str = 'Moving up';
            elseif label == C.G_D
                type_str = 'Moving down';
            elseif label == C.G_NONE
                type_str = 'Other';
            else
                type_str = '(Unrecognized gesture type)';
            end
        end
        
        function type_str = get_wearing_wrist_name(label)
            if label == C.W_L
                type_str = 'Left wrist';
            elseif label == C.W_L
                type_str = 'Right wrist';
            elseif label == C.W_R
                type_str = '(unrecognized wrist)';
            end
        end
        
        function label = get_gesture_movement_label(str)
            if strcmp(str, 'L') == 1
                label = C.G_L;
            elseif strcmp(str, 'R') == 1
                label = C.G_R;
            elseif strcmp(str, 'U') == 1
                label = C.G_U;
            elseif strcmp(str, 'D') == 1
                label = C.G_D;
            else
                error('Unrecognized gesture symbol');
            end
        end
        
        function label = get_wearing_wrist_label(str)
            if strcmp(str, 'L') == 1
                label = C.W_L;
            elseif strcmp(str, 'R') == 1
                label = C.W_R;
            else
                error('Unrecognized wearing wrist symbol');
            end
        end
    end
end